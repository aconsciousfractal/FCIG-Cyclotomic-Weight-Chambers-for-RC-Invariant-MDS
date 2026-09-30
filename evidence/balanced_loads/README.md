# Balanced stable-PCR loads: selected DNA8 certificate

Run from the repository root with Python 3 and no additional dependencies:

```text
python -B -X utf8 scripts/verify_balanced_loads.py
python -B -O -X utf8 scripts/verify_balanced_loads.py
```

The replay reconstructs the DNA8 carrier and the old pool of 1804 mixed walks,
1031 rows, and 622 stable coordinates. It checks the 8102-arc selected core,
exact longest path 106, and a signed rational lift of mass 128, counting each
stable arc with its multiplicity. Forty returns on five stable PCRs certify
stable load at least 133 and full load at least 8235 within this core's real
E0-gradient class. The manuscript proves the transfer to all such gradients;
the replay does not sample potentials. It establishes neither an endpoint nor
an attained class minimum nor a global lower bound.

The seven imported payloads are unchanged; their hashes and source roles are
in PROVENANCE.json. The exact kernel is a shared dependency, not an independent
second graph implementation. The wrapper copies the capsule to an isolated
temporary directory and does not modify shipped inputs. MANIFEST_SHA256.txt
binds all capsule files except itself.
