# Imports
import os
from tqdm.auto import tqdm
import shlex
import csv
from multiprocessing import Pool
import subprocess

OUTDIR = "jobs/2026-01-17.failed"

def unpigz(file):
    command = shlex.split(f"pigz -dkc {file}")
    decompressed = subprocess.Popen(command, stdout=subprocess.PIPE)
    chunk_size = 1024**3 # 1 GiB

    basename= os.path.basename(file).rstrip(".gz") 
    with open(os.path.join(OUTDIR, basename), "wb") as f:
        while True:
            chunk = decompressed.stdout.read(chunk_size) # type: ignore
            if not chunk:
                break
            f.write(chunk)

def main():
    with open(os.path.join(OUTDIR, "metadata.tsv"), "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        print(reader.fieldnames)
        files = []
        for row in reader:
            files.append(row["File1"])
            files.append(row["File2"])

            # f1_gb = os.path.getsize(row["File1"]) / 1000**3
            # f2_gb = os.path.getsize(row["File2"]) / 1000**3
            # print(f"File sizes for Barcode: {row['Barcode']}")
            # print(f"File1 size: {f1_gb} GB")
            # print(f"File2 size: {f2_gb} GB")

        # files = files[:1] # DEBUG

    with Pool(processes=16) as pool:
        for _ in tqdm(pool.imap_unordered(unpigz, files), total=len(files)):
            continue

if __name__ == "__main__":
    main()

# File sizes for Barcode: EB100091
# File1 size: 36.606041562 GB
# File2 size: 37.995361457 GB
# File sizes for Barcode: EB100207
# File1 size: 34.515667817 GB
# File2 size: 36.444314403 GB
# File sizes for Barcode: EB100303
# File1 size: 36.413718935 GB
# File2 size: 38.705112764 GB
# File sizes for Barcode: EB100342
# File1 size: 36.558958222 GB
# File2 size: 38.62846827 GB
# File sizes for Barcode: EB100516
# File1 size: 37.5450933 GB
# File2 size: 39.307683725 GB
# File sizes for Barcode: EB100523
# File1 size: 36.909590496 GB
# File2 size: 38.542559143 GB
# File sizes for Barcode: EB100534
# File1 size: 38.01016198 GB
# File2 size: 39.697892333 GB
