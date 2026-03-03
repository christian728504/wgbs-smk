import os
import sys
import csv
import shutil


def main():

    root_dir = "/data/zusers/ramirezc/wgbs-smk"
    bam_path_fmt = "results/merge/{acc}/{acc}.bam" 
    new_bam_path_fmt = "results/mapping/{acc}/{acc}.bam"

    with open(os.path.join(root_dir, "jobs/2026-02-14/metadata.tsv"), "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        try:
            match sys.argv[1]:
                case "move":
                    for row in reader:
                        acc = row["Barcode"]
                        bam_path = os.path.join(root_dir, bam_path_fmt.format(acc=acc))
                        new_bam_path = os.path.join(root_dir, new_bam_path_fmt.format(acc=acc))

                        if os.path.exists(new_bam_path):
                            print(f"File already exists: {new_bam_path}")
                            continue

                        shutil.move(bam_path, new_bam_path)
                case "copy":
                    for row in reader:
                        acc = row["Barcode"]
                        bam_path = os.path.join(root_dir, bam_path_fmt.format(acc=acc))
                        new_bam_path = os.path.join(root_dir, new_bam_path_fmt.format(acc=acc))

                        if os.path.exists(new_bam_path):
                            print(f"File already exists: {new_bam_path}")
                            continue

                        shutil.copy(bam_path, new_bam_path)
                case "dryrun":
                    for row in reader:
                        acc = row["Barcode"]
                        bam_path = os.path.join(root_dir, bam_path_fmt.format(acc=acc))
                        new_bam_path = os.path.join(root_dir, new_bam_path_fmt.format(acc=acc))
                        
                        if os.path.exists(new_bam_path):
                            print(f"File already exists: {new_bam_path}")
                            continue

                        print(f"mv {bam_path} {new_bam_path}")
                case _:
                    raise ValueError(f"Invalid argument: {sys.argv[1]}")
        except IndexError:
            raise ValueError("Missing argument")

if __name__ == "__main__":
    main()
