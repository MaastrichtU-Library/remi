# REMI 2.1 by Maastricht University Library

## Licence

The REMI source code is licensed under the
[GNU Affero General Public License v3.0](LICENSE).

Documentation and original visual materials are licensed under
[CC BY-NC-SA 4.0](LICENSE-CONTENT.md).

Maastricht University names, logos, trademarks, and other third-party materials are excluded from these licences.

REMI 2.1 (Related Expression Machine Identification) is a web application that helps identify related expressions for PubMed searches.This tool helps researchers explore keyword variants and proximity operators in PubMed searches efficiently.
This short guide explains how to run REMI on your own computer. No programming experience is required.

## What you need

- A Windows, macOS, or Linux computer
- An internet connection, because REMI sends searches to PubMed/NCBI
- Python 3.10 or 3.11
- A free NCBI account

The project is configured for Python 3.10. Python 3.11 is recommended for a new local installation. Do not use a Python version older than 3.10.

## 1. Download REMI from GitHub

1. Open the REMI repository on GitHub.
2. Select the green **Code** button.
3. Select **Download ZIP**.
4. Open your Downloads folder.
5. Right-click the downloaded ZIP file and choose **Extract All** on Windows, or double-click it on macOS.
6. Open the extracted folder.

You are in the correct folder when you can see at least these files:

```text
main.py
config.py
config.ini
requirements.txt
```

Keep all extracted files and folders together. REMI will not work if only `main.py` is copied elsewhere.

## 2. Install Python

Download Python from the [official Python website](https://www.python.org/downloads/).

Python 3.11 is recommended for this version of REMI.

### Windows

1. Run the Python installer.
2. On the first installer screen, select **Add python.exe to PATH**.
3. Complete the installation.

Open **Command Prompt** and check the installed version:

```bat
py -3.11 --version
```

The result should start with `Python 3.11`.

### macOS or Linux

Open **Terminal** and check the installed version:

```bash
python3 --version
```

The result must be Python 3.10 or 3.11. If Python is missing or too old, install Python 3.11 before continuing.

## 3. Open a command window in the REMI folder

All commands in the rest of this guide must be run from the folder containing `main.py`.

### Windows

1. Open the extracted REMI folder in File Explorer.
2. Select the address bar at the top of File Explorer.
3. Type `cmd` and press Enter.

A Command Prompt window opens in the correct folder.

### macOS or Linux

Open Terminal, type `cd` followed by a space, drag the extracted REMI folder into the Terminal window, and press Enter. For example:

```bash
cd "/path/to/the/REMI-folder"
```

## 4. Create a virtual environment

A virtual environment keeps REMI's Python packages separate from other programs on your computer. It only needs to be created once.

### Windows — Command Prompt

```bat
py -3.11 -m venv .venv
.venv\Scripts\activate.bat
```

### macOS or Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
```

When activation succeeds, the command line normally begins with `(.venv)`.

## 5. Install the required packages

With the virtual environment active, run these two commands:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The installation may take a few minutes. Warnings are usually harmless, but any error must be resolved before continuing.

## 6. Create the Flask secret key

The `SECRET_KEY` protects REMI's browser session cookies. Generate a new random key by running:

```bash
python -c "import secrets; print(secrets.token_hex(32))"
```

The command prints a long line of letters and numbers. Copy that complete line. Each REMI 2.1 installation should have its own secret key.

## 7. Obtain an NCBI API key

The `NCBI_API_KEY` identifies REMI's requests to PubMed and provides the appropriate NCBI request allowance.

1. Sign in to your [NCBI account](https://account.ncbi.nlm.nih.gov/).
2. Open **Account settings**.
3. Find **API Key Management**.
4. Select **Create an API Key**.
5. Copy the generated key.

NCBI explains API-key use in its [official E-utilities documentation](https://www.ncbi.nlm.nih.gov/books/NBK25497/).

## 8. Fill in `config.ini`

The configuration file is named `config.ini`. It must remain in the same folder as `main.py`.

Open `config.ini` in a plain-text editor such as Notepad. Its contents should look like this:

```ini
[REMI]

# Random Flask key used to securely sign browser session cookies.
# Generate it with: python -c "import secrets; print(secrets.token_hex(32))"
SECRET_KEY = paste_your_new_secret_key_here

# Personal NCBI API key used for PubMed requests.
# Create it in the API Key Management section of your NCBI account settings.
NCBI_API_KEY = paste_your_ncbi_api_key_here
```

Replace both placeholder values with the keys you created. Do not add quotation marks around the keys.

Save the file and make sure your editor has not renamed it to `config.ini.txt`.

> **Keep these keys private.** Do not publish your completed `config.ini`, commit it to GitHub, send it by email, or include it in a public ZIP file. A public repository should contain only a placeholder/example configuration file.

## 9. Start REMI

Make sure the virtual environment is still active, then run:

```bash
python -m flask --app main run --host 127.0.0.1 --port 8080
```

When REMI has started, open this address in a web browser:

[http://127.0.0.1:8080](http://127.0.0.1:8080)

Keep the command window open while using REMI. To stop the app, return to the command window and press **Ctrl+C**.

This command makes REMI available only on the computer on which it is running. The built-in Flask server is intended for local use and testing, not for publishing REMI on the public internet.

## 10. Start REMI again later

The next time you want to use REMI:

1. Open a command window in the REMI folder as described in step 3.
2. Activate the existing virtual environment.

On Windows:

```bat
.venv\Scripts\activate.bat
```

On macOS or Linux:

```bash
source .venv/bin/activate
```

3. Start the app:

```bash
python -m flask --app main run --host 127.0.0.1 --port 8080
```

4. Open [http://127.0.0.1:8080](http://127.0.0.1:8080) in your browser.

You do not need to recreate the virtual environment or reinstall the packages each time.

## Troubleshooting

### “Python was not found” or “py is not recognized”

Reinstall Python 3.11 and make sure **Add python.exe to PATH** is selected in the Windows installer. Close and reopen Command Prompt afterward.
With a mac any python command may be "python3" instead of "python".

### PowerShell says that running scripts is disabled

The Windows instructions above use Command Prompt and avoid this issue. Open the REMI folder in File Explorer, type `cmd` in the address bar, and use the Command Prompt commands instead.

### “No module named …”

The virtual environment is probably inactive, or the packages were not installed. Activate `.venv` and run:

```bash
python -m pip install -r requirements.txt
```

### An error mentions `config.ini`, `SECRET_KEY`, or `NCBI_API_KEY`

Check that:

- The file is named exactly `config.ini`.
- It is in the same folder as `main.py`.
- The first line is `[REMI]`.
- Both key values are present.
- You started REMI from the folder containing `main.py`.

### Port 8080 is already in use

Another program is using the same port. Start REMI on port 8081 instead:

```bash
python -m flask --app main run --host 127.0.0.1 --port 8081
```

Then open [http://127.0.0.1:8081](http://127.0.0.1:8081).

### PubMed requests fail

Check your internet connection and confirm that the NCBI key in `config.ini` is complete and belongs to your NCBI account. If the key was ever published or shared, revoke it in your NCBI account and create a new one.

## Important security notes

- Never share the real values from `config.ini`.
- Every person installing REMI should generate their own `SECRET_KEY` and obtain their own `NCBI_API_KEY`.
- If either key is exposed, replace it immediately.
- Do not use the local Flask server to host a public website.
- Stop REMI with **Ctrl+C** when you have finished using it.
