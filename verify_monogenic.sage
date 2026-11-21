# Verifies generators of monogenic fields
# If all asserts pass, then everything is correct

from lmf import db    # or with lmfdb-lite:  from lmf import db

R.<x> = QQ[]

# Read in data from file
print("Reading in data from file...")
lines = open("monogenic_coeffs_v1.txt").read().split("\n")
rows = [line.split("|") for line in lines[3:] if line]

print("Verifying given generators...")
for i in range(0, len(rows), 10000):   # ask the LMFDB for 10000 fields at a time
    batch = rows[i:i+10000]
    lmfdb = {rec["label"]: rec for rec in db.nf_fields.search({"label": {"$in": [label for label, _ in batch]}}, ["label", "coeffs", "disc_sign", "disc_abs"])}
    for label, coeffs in batch:
        K = lmfdb[label]                
        n, r, D, _ = [ZZ(c) for c in label.split(".")]

        f = R([ZZ(int(c)) for c in coeffs[1:-1].split(",")])   # our claimed monogenic generator
        g = R([ZZ(int(c)) for c in K["coeffs"]])               # the LMFDB polynomial of K

        # 1. Check f (and g) is monic and irreducible of degree n
        assert f.degree()==n and f.is_monic() and f.is_irreducible(), label
        assert g.degree()==n and g.is_monic() and g.is_irreducible(), label

        # 2. Check Q(alpha) = K:
        L.<a> = NumberField(f)
        M.<b> = NumberField(g)
        assert L.is_isomorphic(M), label

        # 3. Check disc(f) = disc(K), from the LMFDB and from the label.
        disc_K = ZZ(int(K["disc_sign"])) * ZZ(int(K["disc_abs"]))
        assert f.discriminant()==disc_K, label

    print(i+len(batch), "fields checked, last label:", label, flush=True)

print("All", len(rows), "given generators are correct! :)")
