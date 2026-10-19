# config.py
import configparser
import os

config = configparser.ConfigParser()
# config.ini sits in the same folder as this file
config.read(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.ini'))

def _is_placeholder(value):
    v = (value or '').strip().lower()
    return (not v or v == 'fallback-secret' or v.startswith('[')
            or 'paste' in v or 'fill in' in v or 'your_' in v)

SECRET_KEY = config.get('REMI', 'SECRET_KEY', fallback='')
if _is_placeholder(SECRET_KEY):
    raise RuntimeError(
        "SECRET_KEY in config.ini is missing or still a placeholder. Generate one with\n"
        '    python -c "import secrets; print(secrets.token_hex(32))"\n'
        "and paste it into config.ini (see README, steps 6 and 8)."
    )

# Without an API key REMI still works, but NCBI then allows only 3 requests per second.
NCBI_API_KEY = config.get('REMI', 'NCBI_API_KEY', fallback='').strip()
if _is_placeholder(NCBI_API_KEY):
    NCBI_API_KEY = ''

# Optional contact address that NCBI can use if there is a problem with REMI's requests.
CONTACT_EMAIL = config.get('REMI', 'CONTACT_EMAIL', fallback='').strip()
if _is_placeholder(CONTACT_EMAIL):
    CONTACT_EMAIL = ''
