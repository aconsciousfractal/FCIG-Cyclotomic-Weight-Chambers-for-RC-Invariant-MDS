# Reproduce version 3.2.0

Requires Python 3.12 or newer and requirements.txt. All scientific data are local.

    python -m pip install -r requirements.txt
    python -X utf8 -B scripts/verify.py --profile manifest
    python -X utf8 -B scripts/verify.py --profile core
    python -X utf8 -B scripts/verify_extended.py --optimized
    python -X utf8 -B scripts/verify_advances.py --optimized
    python -X utf8 -B scripts/verify_balanced_loads.py --optimized

Core includes the normal extension and advances replays, plus both modes of the balanced-load replay. The latter also has five semantic mutation controls in tests/test_balanced_loads.py and runs alone without third-party Python packages. The new runner forwards -O to every child explicitly; assertion-based imported tests have been made explicit exceptions. Use --logs OUT on verify_advances.py to retain detailed per-check output, and --output OUT.json for a summary. Temporary copies preserve committed receipts. UTF-8/LF and POSIX manifest keys make comparisons portable.

On a clean complete Git checkout, full also verifies isolated materialization and the one-root history. Normal revision commits are allowed.

Install Tectonic 0.16.9 for the PDF, then run:

    python -X utf8 -B scripts/build_paper.py --verify-rebuild

TECTONIC may specify the executable. The fixed SOURCE_DATE_EPOCH is recorded. Initial TeX bundle caching may need internet; scientific checks do not.

The extended README explains finite domains and independence. No timeout is interpreted as a result. No VPS or remote execution is configured.
