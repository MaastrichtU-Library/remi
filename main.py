import xml.etree.ElementTree as ET
import requests
import time
import json
import re
import os
from collections import Counter
from flask import Flask, request, render_template, jsonify, Response, stream_with_context, session
from metapub import PubMedFetcher
import uuid
from datetime import datetime, timedelta
import shutil
import threading
from functools import wraps
import logging
from collections import defaultdict
import config  # this loads config.py, which handles the two secrets in config.txt

_cleanup_lock = threading.Lock()
_last_cleanup_ts = 0

ALLOWED_RE = re.compile(r'^[A-Za-z0-9*\s]+$')

def norm_ws(s: str) -> str:
    # collapse whitespace only, do NOT touch hyphens
    return re.sub(r'\s+', ' ', s.strip().lower())

def validate_query_text(raw: str) -> str:
    if raw is None:
        raise ValueError("Empty input.")

    text = raw.strip()

    # 1) Only ASCII alnum, whitespace, and *
    if not ALLOWED_RE.match(text):
        raise ValueError("Only A–Z, 0–9, spaces/newlines, and * are allowed. No other characters.")

    # 2) Validate each token (split on whitespace/newlines)
    for token in re.split(r'\s+', text):
        if not token:
            continue

        # PubMed requires >=4 characters before the first wildcard
        if '*' in token:
            prefix = token.split('*', 1)[0]
            if len(prefix) < 4:
                raise ValueError(f"Wildcard terms must have at least 4 characters before the first * (problem: '{token}').")

    return text

CLEANUP_INTERVAL_SECONDS = 15 * 60  # run cleanup at most every 15 minutes

# set-up the secret keys. Make sure you have set them in the config.txt
app = Flask(__name__)
app.config['SECRET_KEY'] = config.SECRET_KEY
os.environ['NCBI_API_KEY'] = config.NCBI_API_KEY

# Global query throttle
QUERY_LOCK = threading.Lock()
GLOBAL_THROTTLE_DELAY = 0.12  # max ~8 queries/sec

# Per-session pacing
MIN_SESSION_DELAY = 1.0  # seconds between same-session queries

class PartialResultsError(Exception):
    def __init__(self, message, partial_results):
        super().__init__(message)
        self.partial_results = partial_results
        
def maybe_run_cleanup():
    """Run cleanup infrequently to avoid disk-walk latency on every request."""
    global _last_cleanup_ts
    now = time.time()

    # Fast-path: avoid lock for the common case
    if now - _last_cleanup_ts < CLEANUP_INTERVAL_SECONDS:
        return

    # Only one request should actually do the work
    with _cleanup_lock:
        now = time.time()
        if now - _last_cleanup_ts < CLEANUP_INTERVAL_SECONDS:
            return

        cleanup_old_sessions()
        cleanup_old_cache()

        _last_cleanup_ts = now

@app.context_processor
def inject_session_id():
    return {'session_id': session.get('id', 'N/A')}

def touch_session(sid: str) -> None:
    """Mark a session as active (for cleanup/expiry)."""
    if not sid:
        return
    user_progress.setdefault(sid, {})
    user_progress[sid]['last_seen'] = datetime.utcnow()

alt_cache_dir = '/tmp/eutils_cache'
user_progress = {}  # session_id: { search_data..., 'last_seen': datetime }
USER_TIMEOUT = timedelta(minutes=30)
CACHE_DIR = '/tmp/eutils_cache'
CACHE_MAX_AGE_DAYS = 30

@app.before_request
def ensure_session_id():
    if 'id' not in session:
        session['id'] = str(uuid.uuid4())

    sid = session['id']

    # If this session has progress tracked, keep it alive
    if sid in user_progress:
        touch_session(sid)

    # print(f"[DEBUG] Active session ID: {sid}")
    maybe_run_cleanup()

def cleanup_old_sessions():
    now = datetime.utcnow()

    expired_ids = []
    for sid, entry in list(user_progress.items()):
        last_seen = entry.get('last_seen')
        if not isinstance(last_seen, datetime):
            # If last_seen is missing or malformed, expire it
            expired_ids.append(sid)
            continue

        if now - last_seen > USER_TIMEOUT:
            expired_ids.append(sid)

    for sid in expired_ids:
        # print(f"[DEBUG] Cleaning up expired session: {sid}")
        del user_progress[sid]

def cleanup_old_cache():
    now = time.time()
    max_age = CACHE_MAX_AGE_DAYS * 86400  # seconds in 30 days

    if not os.path.exists(CACHE_DIR):
        return

    for root, dirs, files in os.walk(CACHE_DIR):
        for name in files:
            filepath = os.path.join(root, name)
            try:
                if os.path.isfile(filepath):
                    file_age = now - os.path.getmtime(filepath)
                    if file_age > max_age:
                        os.remove(filepath)
                        #print(f"[DEBUG] Deleted old cache file: {filepath}")
            except Exception as e:
                print(f"[DEBUG] Error deleting cache file {filepath}: {e}")

def find_precise_variants(term1_list, term2_list, date_from=None, date_to=None, exclude_hyphens=False):
    print("[DEBUG] ENTERED find_precise_variants()")
    all_variants = {'term1': {}, 'term2': {}}

    total = len(term1_list) * len(term2_list)
    completed = 0

    for t1 in term1_list:
        for t2 in term2_list:
            print(f"[DEBUG] Running precise comparison: {t1} AND {t2}")
            
            # Run the variant finder (may have many internal iterations)
            result = find_combined_variants(t1, t2, date_from, date_to, exclude_hyphens)

            # Merge results into all_variants
            for key in ['term1', 'term2']:
                for base, variants in result[key].items():
                    if base not in all_variants[key]:
                        all_variants[key][base] = variants
                    else:
                        all_variants[key][base].update(variants)

            # Now update precise_progress — once per comparison (not per iteration)
            completed += 1
            if session['id'] not in user_progress:
                user_progress[session['id']] = {}

            user_progress[session['id']]['precise_progress'] = {
                'current': completed,
                'total': total,
                'last_updated': datetime.utcnow().isoformat()
            }

    return all_variants

# --- HY-LOGIC PATCH for hyphen-exclusion behavior ---
def handle_hyphen_variant(match_lower, raw_seen_prefixes, seen_variants, found_variants):
    prefix = match_lower.split('-')[0]
    seen_variants.add(prefix)
    if prefix in raw_seen_prefixes:
        found_variants[prefix] += 1
    else:
        raw_seen_prefixes.add(prefix)
        found_variants[prefix] += 1
                
progress = {}

# ---------- UTILITIES ----------
def rate_limited_request(fn):
    @wraps(fn)
    def wrapper(*args, **kwargs):
        # Session pacing
        last = session.get('last_pubmed_call', 0)
        now = time.time()
        since = now - last
        if since < MIN_SESSION_DELAY:
            time.sleep(MIN_SESSION_DELAY - since)
        session['last_pubmed_call'] = time.time()

        # Global throttle
        with QUERY_LOCK:
            time.sleep(GLOBAL_THROTTLE_DELAY)
            return fn(*args, **kwargs)
    return wrapper

def format_terms(terms):
    # Match either quoted phrases or single terms
    pattern = r'"[^"]+"|\S+'
    parts = re.findall(pattern, terms)
    clean_parts = []
    buffer = []

    for part in parts:
        if part.lower() == 'or' and buffer:
            clean_parts.append(' '.join(buffer))
            buffer = []
        else:
            buffer.append(part)

    if buffer:
        clean_parts.append(' '.join(buffer))

    formatted = []
    for term in clean_parts:
        stripped = term.strip()
        if stripped.startswith('"') and stripped.endswith('"'):
            # Already quoted exact phrase
            formatted.append(f'{stripped}[tiab]')
        else:
            # Unquoted term or phrase, quote it first
            if ' ' in stripped:
                stripped = f'"{stripped}"'
            formatted.append(f'{stripped}[tiab]')

    return " OR ".join(formatted)

def format_combined_query(term1, term2):
    return f"({format_terms(term1)}) AND ({format_terms(term2)})"

def update_combined_query(term1, term2, term1_exclusions, term2_exclusions):
    def format_exclusions(exclusions):
        if not exclusions:
            return ""
        return " NOT (" + " OR ".join([f'\"{e}\"[tiab]' for e in sorted(exclusions)]) + ")"

    term1_clause = f"({term1}){format_exclusions(term1_exclusions)}"
    term2_clause = f"({term2}){format_exclusions(term2_exclusions)}"

    return f"({term1_clause}) AND ({term2_clause})"

@rate_limited_request
def get_pmids(query, retmax=1000, date_from=None, date_to=None):
    from metapub import PubMedFetcher
    fetch = PubMedFetcher(cachedir=alt_cache_dir)

    date_range = ""
    if date_from and date_to:
        date_range = f' AND ("{date_from}"[PDAT] : "{date_to}"[PDAT])'
    elif date_from:
        date_range = f' AND ("{date_from}"[PDAT] : "3000"[PDAT])'
    elif date_to:
        date_range = f' AND ("1800"[PDAT] : "{date_to}"[PDAT])'

    full_query = query + date_range

    last_exception = None

    for attempt in range(3):
        try:
            return fetch.pmids_for_query(query=full_query, retmax=retmax)
        except Exception as e:
            last_exception = e
            delay = 2 ** attempt
            print(f"[DEBUG] get_pmids() failed (attempt {attempt + 1}/3): {e}")
            time.sleep(delay)

    final_error = f"get_pmids() failed after 3 attempts. Last error: {last_exception}"
    print(f"[DEBUG] {final_error}")
    logging.error(final_error)
    raise Exception(final_error)

@rate_limited_request
def extract_title_abstract(pmid_chunk):
    import xml.etree.ElementTree as ET
    import requests

    pmid_str = ','.join(pmid_chunk)
    efetch = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id={pmid_str}"
    headers = {"User-Agent": "REMI/2.0 (rinus.verdonschot@maastrichtuniversity.nl)"}

    last_exception = None

    for attempt in range(3):
        try:
            response = requests.get(efetch, headers=headers, timeout=10)
            response.raise_for_status()

            root = ET.fromstring(response.content)
            articles = []

            for article in root.findall("PubmedArticle"):
                pmid = article.find("MedlineCitation/PMID").text
                title = article.find("MedlineCitation/Article/ArticleTitle").text or ""
                abstract_section = article.find("MedlineCitation/Article/Abstract")
                abstract_text = ""

                if abstract_section is not None:
                    pieces = []
                    for elem in abstract_section.findall("AbstractText"):
                        part = elem.text or ""
                        label = elem.attrib.get("Label")
                        if label:
                            pieces.append(f"{label}: {part}")
                        else:
                            pieces.append(part)
                    abstract_text = ' '.join(pieces).strip()
                else:
                    abstract_text = "No abstract"

                keywords = [kw.text for kw in article.findall(".//Keyword") if kw.text]
                full_text = f"{title} {abstract_text} {' '.join(keywords)}".lower()

                articles.append({"pmid": pmid, "full_text": full_text})

            return articles

        except requests.RequestException as e:
            last_exception = e
            delay = 2 ** attempt
            msg = f"PubMed fetch failed (try {attempt + 1}/3): {e}"
            print(f"[DEBUG] {msg}")
            logging.warning(msg)
            time.sleep(delay)

        except ET.ParseError as e:
            msg = f"XML parsing failed: {e}"
            print(f"[DEBUG] {msg}")
            logging.warning(msg)
            raise Exception("PubMed response could not be parsed. Please try again later.")

        except Exception as e:
            msg = f"Unexpected error in extract_title_abstract(): {e}"
            print(f"[DEBUG] {msg}")
            logging.error(msg)
            raise Exception("Unexpected error occurred while processing PubMed response.")

    if last_exception:
        final_error = f"PubMed fetch failed after 3 attempts: {last_exception}"
    else:
        final_error = "PubMed fetch failed after 3 attempts due to unknown error."

    print(f"[DEBUG] {final_error}")
    logging.error(final_error)
    raise Exception(final_error)

def chunk_pmids(pmids, chunk_size=200):
    for i in range(0, len(pmids), chunk_size):
        yield pmids[i:i + chunk_size]

def format_variants_for_display(variants_dict):
    formatted = {}
    for term, counter in variants_dict.items():
        formatted[term] = []
        for word, count in counter.most_common():
            if word.startswith('__H__'):
                clean = word.replace('__H__', '')
                display = f"{clean} [+hyphen]"
            else:
                display = word
            count_display = "non-wildcard" if count == 9999 else count
            formatted[term].append({"word": display, "count": count_display})
    return formatted

# ---------- KEYWORD VARIANT LOGIC ----------

def find_keyword_variants(term, date_from=None, date_to=None, exclude_hyphens=False):
    term_clean = term.strip().lower()

    # Build a regex that treats '*' as "0+ word chars"
    if '*' in term_clean:
        # PubMed wildcard matches 0 or more characters within a term.
        # We'll approximate that with \w* (letters/digits/_), and allow hyphens too.
        # Escape everything else, then replace escaped \* with a pattern.
        pattern_core = re.escape(term_clean).replace(r'\*', r'[\w-]*')
        pattern = rf"\b{pattern_core}\b"
    else:
        base_term = term_clean
        pattern = rf"\b{re.escape(base_term)}[\w-]*\b"

    found_variants = Counter()
    seen_variants = set()
    raw_seen_prefixes = set()

    query = f"{term.strip()}[tiab]"

    iteration = 0
    try:
        while iteration < 10:
            pmids = get_pmids(query, date_from=date_from, date_to=date_to)
            if not pmids:
                break

            new_found = False
            for pmid_chunk in chunk_pmids(pmids):
                articles = extract_title_abstract(pmid_chunk)
                for article in articles:
                    text = article.get('full_text', '').lower()
                    matches = re.findall(pattern, text)

                    for match in matches:
                        match_lower = match.lower()
                        is_hyphenated = '-' in match_lower

                        if not is_hyphenated:
                            if match_lower not in seen_variants:
                                seen_variants.add(match_lower)
                                found_variants[match_lower] += 1
                                raw_seen_prefixes.add(match_lower)
                                new_found = True
                            else:
                                found_variants[match_lower] += 1

                        elif exclude_hyphens:
                            handle_hyphen_variant(match_lower, raw_seen_prefixes, seen_variants, found_variants)
                            new_found = True
                        else:
                            if match_lower not in seen_variants:
                                seen_variants.add(match_lower)
                                found_variants[match_lower] += 1
                                new_found = True
                            else:
                                found_variants[match_lower] += 1

            if not new_found:
                break

            clean_exclusions = {v.replace('__H__', '') for v in seen_variants}
            exclusions = " NOT (" + " OR ".join([f'\"{v}\"[tiab]' for v in sorted(clean_exclusions)]) + ")"
            query = f"{term.strip()}[tiab]{exclusions}"
            print(f"[DEBUG] Iteration {iteration + 1} - PubMed query:\n{query}\n")
            iteration += 1

        print("[DEBUG] Found variants:", found_variants)
        return [{"word": w.replace('__H__', ''), "count": c} for w, c in found_variants.most_common()]

    except Exception as e:
        partials = [{"word": w.replace('__H__', ''), "count": c} for w, c in found_variants.most_common()]
        raise PartialResultsError(str(e), partials)


# ---------- COMBINED VARIANT LOGIC ----------
def find_combined_variants(term1_query, term2_query, date_from=None, date_to=None, exclude_hyphens=False):
    variants = {'term1': {}, 'term2': {}}
    discovered_variants = {'term1': set(), 'term2': set()}
    raw_seen_prefixes = {'term1': set(), 'term2': set()}
    truncated_terms = {'term1': [], 'term2': []}
    wildcard_seeds = {'term1': set(), 'term2': set()}  # ⬅️ key fix
    
    literal_seeds = {'term1': set(), 'term2': set()}  # exact (non-wildcard) user terms
    
    seed_to_base = {'term1': {}, 'term2': {}}

    for key in ['term1', 'term2']:
        query_string = term1_query if key == 'term1' else term2_query
        for t in re.split(r'\s+or\s+', query_string, flags=re.I):
            t = t.strip()
            if not t:
                continue
            if t.endswith('*'):
                seed = t.lower()  # keep full wildcard seed e.g. "tumo*r*"
                base = t.rstrip('*').lower()  # key name in variants dict (we can keep this)
                variants[key][base] = Counter()
                truncated_terms[key].append(seed)          # store seed, not base
                wildcard_seeds[key].add(seed)
                seed_to_base[key][seed] = base

            else:
                word = t.lower().strip()
                word_norm = norm_ws(word)   # normalize spaces/hyphens same way as matches                            
                variants[key][word_norm] = Counter({word_norm: 9999})

                # Track literal seeds separately so they won't be excluded later
                literal_seeds[key].add(word_norm)

                # IMPORTANT: do NOT add literal seeds to discovered_variants,
                # because discovered_variants are used to build NOT-exclusions.

    def format_truncated_only(query_string):
        terms = [t.strip() for t in re.split(r'\s+or\s+', query_string, flags=re.I) if t.strip()]
        return " OR ".join([f"{t}[tiab]" for t in terms])

    def build_exclusion_clause(all_discovered):
        if not all_discovered:
            return ""

        clean_terms = []
        for w in sorted(all_discovered):
            clean = w.replace('__H__', '')
            if clean.endswith('*'):
                continue
            if clean.startswith('"') and clean.endswith('"'):
                clean = clean[1:-1]
            clean_terms.append(f'"{clean}"[tiab]')

        if not clean_terms:
            return ""
        return " NOT (" + " OR ".join(clean_terms) + ")"

    term1_base = format_truncated_only(term1_query)
    term2_base = format_truncated_only(term2_query)
    term1_exclusions = set()
    term2_exclusions = set()

    iteration = 0
    while iteration < 21:
        full_query = (
            f"(({term1_base}){build_exclusion_clause(term1_exclusions)}) "
            f"AND "
            f"(({term2_base}){build_exclusion_clause(term2_exclusions)})"
        )
        print(f"[DEBUG] Iteration {iteration + 1} - Combined PubMed query:\n{full_query}\n")

        pmids = get_pmids(full_query, date_from=date_from, date_to=date_to)
        if not pmids:
            break

        new_found = False
        for pmid_chunk in chunk_pmids(pmids):
            articles = extract_title_abstract(pmid_chunk)
            for article in articles:
                text = article['full_text']
                for key, term_list in [('term1', truncated_terms['term1']), ('term2', truncated_terms['term2'])]:
                    for seed in term_list:
                        base = seed_to_base[key].get(seed, seed.rstrip('*'))
                        base_norm = seed.replace('*', '').lower()

                        if ' ' in seed:
                            parts = seed.split()
                            if len(parts) == 2:
                                prefix, wildcard = parts
                                pattern = rf"\b{re.escape(prefix)}[-\s]+{re.escape(wildcard)}[\w-]*\b"
                                matches = re.findall(pattern, text)
                            else:
                                matches = []
                        else:
                            pattern_core = re.escape(seed).replace(r'\*', r'[\w-]*')
                            pattern = rf"\b{pattern_core}\b"
                            matches = re.findall(pattern, text)

                        for match in matches:
                            match_lower = match.lower()
                            match_norm  = norm_ws(match_lower)

                            is_hyphenated = '-' in match_lower

                            # If exclude_hyphens is ON: ALWAYS count the prefix, NEVER the hyphenated form
                            if exclude_hyphens and is_hyphenated:
                                prefix = match_lower.split('-', 1)[0]   # or match_norm.split('-',1)[0]
                                prefix = norm_ws(prefix)

                                if prefix not in discovered_variants[key]:
                                    discovered_variants[key].add(prefix)
                                    new_found = True

                                variants[key][base][prefix] += 1
                                continue

                            # ---- normal behavior (keep hyphenated variants) ----
                            variant_key = match_norm
                            if is_hyphenated and not exclude_hyphens:
                                variant_key = "__H__" + match_norm

                            if variant_key not in discovered_variants[key]:
                                discovered_variants[key].add(variant_key)
                                new_found = True

                            variants[key][base][variant_key] += 1

        if not new_found:
            break

        term1_exclusions = set(discovered_variants['term1']) - literal_seeds['term1']
        term2_exclusions = set(discovered_variants['term2']) - literal_seeds['term2']

        iteration += 1

    return variants


# ---------- ROUTES ----------

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/keyword_search')
def keyword_search():
    return render_template('keyword_search.html')

@app.route('/multiple_search')
def multiple_search():
    return render_template('multiple_search.html')

...
@app.route('/search', methods=['POST'])
def search():
    global progress
    term1_query = request.form['term1']
    term2_query = request.form['term2']
    
    #print("[DEBUG] exclude_hyphens raw:", request.form.get('exclude_hyphens'))
    #print("[DEBUG] exclude_hyphens parsed:", exclude_hyphens)

    # enforce allowed characters + multi-wildcard allowed
    try:
        term1_query = validate_query_text(term1_query)
        term2_query = validate_query_text(term2_query)
    except ValueError as e:
        return jsonify({"status": "error", "message": str(e)}), 400

    date_from = request.form.get('date_from')
    date_to = request.form.get('date_to')
    exclude_hyphens = str(request.form.get('exclude_hyphens', '')).lower() in ('true', 'on', '1', 'yes')
    #exclude_hyphens = request.form.get('exclude_hyphens') == 'true'
    precise = request.form.get('precise_proximity') == 'true'    
    print(f"[DEBUG] Precise proximity enabled? {precise}")
    
    try:
        term1_list = [t.strip() for t in re.split(r'\s+or\s+', term1_query, flags=re.I) if t.strip()]
        term2_list = [t.strip() for t in re.split(r'\s+or\s+', term2_query, flags=re.I) if t.strip()]

        if precise:
            variants = find_precise_variants(term1_list, term2_list, date_from, date_to, exclude_hyphens)
        else:
            variants = find_combined_variants(term1_query, term2_query, date_from, date_to, exclude_hyphens)

        # Store the search progress in the user-specific dictionary
        user_progress[session['id']] = {
            'term1_query': term1_query,
            'term2_query': term2_query,
            'term1_variants': variants['term1'],
            'term2_variants': variants['term2']
        }
        
        touch_session(session['id'])
        
        return jsonify({"status": "started"})

    except Exception as e:
        import traceback
        user_progress[session['id']] = {
            'term1_query': term1_query,
            'term2_query': term2_query,
            'term1_variants': variants['term1'] if 'variants' in locals() else {},
            'term2_variants': variants['term2'] if 'variants' in locals() else {},
            'error': str(e),
            'trace': traceback.format_exc()
        }
        touch_session(session['id'])
        return jsonify({"status": "error", "message": str(e)}), 500

@app.route('/search_stream', methods=['GET'])
def search_stream():
    global progress

    def generate():
        progress = user_progress.get(session['id'])
        if not progress:
            yield "data: No progress found.\n\n"
            yield "data: DONE\n\n"
            return

        yield f"data: Variants for Term 1 ({progress['term1_query']}): {sum(len(v) for v in progress['term1_variants'].values())}\n\n"
        yield f"data: Variants for Term 2 ({progress['term2_query']}): {sum(len(v) for v in progress['term2_variants'].values())}\n\n"

        results = {
            "term1_variants": format_variants_for_display(progress['term1_variants']),
            "term2_variants": format_variants_for_display(progress['term2_variants'])
        }

        yield f"data: {json.dumps(results)}\n\n"
        yield "data: DONE\n\n"

    return Response(stream_with_context(generate()), mimetype='text/event-stream')

@app.route('/generate_proximity_search', methods=['POST'])
def generate_proximity_search():
    data = request.get_json()
    term1_variants = data.get('term1_variants', [])
    term2_variants = data.get('term2_variants', [])
    proximity_distance = data.get('proximity_distance', 3)

    def clean(word):
        return word.replace(' [+hyphen]', '')

    combinations = [
        f'"{clean(t1)} {clean(t2)}"[tiab:~{proximity_distance}]'
        for t1 in term1_variants for t2 in term2_variants
    ]

    proximity_query = " OR ".join(combinations)
    return jsonify({"proximity_query": proximity_query})

...
@app.route('/get_keyword_variants', methods=['POST'])
def get_keyword_variants():
    data = request.get_json()
    term = data.get('term')
    date_from = data.get('date_from')
    date_to = data.get('date_to')
    exclude_hyphens = data.get('exclude_hyphens', False)
    
    try:
        term = validate_query_text(term)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    try:
        variants = find_keyword_variants(term, date_from=date_from, date_to=date_to, exclude_hyphens=exclude_hyphens)
        return jsonify({"variants": variants})
    except PartialResultsError as e:
        return jsonify({
            "error": str(e),
            "partial_results": e.partial_results
        }), 500
    except Exception as e:
        import traceback
        return jsonify({
            "error": str(e),
            "trace": traceback.format_exc(),
            "partial_results": []
        }), 500

@app.route('/memory_game')
def memory_game():
    return render_template('platform_game.html')

@app.route('/contact')
def contact():
    return render_template('contact.html')

@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html'), 404

@app.route('/precise_progress')
def precise_progress():
    sid = session.get('id')
    progress = user_progress.get(sid, {}).get('precise_progress', {})
    return jsonify(progress)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))