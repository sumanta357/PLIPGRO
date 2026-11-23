import os
import time
import pymol as p
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
#contributor : Sumanta K Sahu 
# =========================
# USER INPUT
# =========================
target_ligand = input(
    "\nEnter ligand residue name (example: LIG): "
).strip().upper()

target_metal = input(
    "\nEnter metal residue name (example: CU1, ZN, FE, MG) [Press Enter if none]: "
).strip().upper()

# =========================
# OUTPUT FILE
# =========================
if target_metal:
    outfile_name = f"{target_ligand}_{target_metal}_md_interactions.csv"
else:
    outfile_name = f"{target_ligand}_md_interactions.csv"

# =========================
# FILES
# =========================
directory_traj = os.getcwd()

traj_file = os.path.join(directory_traj, "md_clean.xtc")
str_file = os.path.join(directory_traj, "md_clean.pdb")

start_time = time.time()

# =========================
# LOAD TRAJECTORY
# =========================
p.cmd.load(str_file)
p.cmd.load_traj(traj_file)

n_frames = p.cmd.count_frames()

# =========================
# FRAME DIRECTORY
# =========================
frames_dir = os.path.join(directory_traj, "frames")
os.makedirs(frames_dir, exist_ok=True)

interaction_counts = defaultdict(int)
frame_interactions = {}

# =========================
# ALLOWED PLIP BINDING SITES
# =========================
allowed_sites = [target_ligand]

if target_metal:
    allowed_sites.append(target_metal)

print("\nAnalyzing binding sites:")
for site in allowed_sites:
    print(f"  {site}")

# =========================
# LOOP OVER FRAMES
# =========================
for i in range(1, n_frames + 1):

    p.cmd.frame(i)

    print(f"\nSaving frame {i}")

    pdb_file = os.path.join(
        frames_dir,
        f"frame_{i:04d}.pdb"
    )

    p.cmd.save(pdb_file, "all")

    frame_id = os.path.splitext(
        os.path.basename(pdb_file)
    )[0]

    print(f"Checking interactions in frame {i}")

    # =========================
    # RUN PLIP
    # =========================
    os.system(
        f'plip -f "{pdb_file}" -x -q'
    )

    report_path = f"{frame_id}_report.xml"

    if not os.path.exists(report_path):
        print(
            f"No PLIP report for frame {i}, skipping"
        )
        continue

    # =========================
    # PARSE XML
    # =========================
    try:
        tree = ET.parse(report_path)
        root = tree.getroot()
    except Exception as e:
        print(f"XML parse error: {e}")
        continue

    this_frame = set()

    binding_sites = root.findall(".//bindingsite")

    # =========================
    # PROCESS BINDING SITES
    # =========================
    for site in binding_sites:

        hetid = site.findtext(
            "./identifiers/hetid"
        )

        if hetid is None:
            continue

        hetid = hetid.strip().upper()

        # -------------------------
        # ONLY PROCESS REQUESTED
        # LIGAND / METAL SITES
        # -------------------------
        if hetid not in allowed_sites:
            continue

        interactions = site.find("interactions")

        if interactions is None:
            continue

        # =========================
        # INTERACTIONS
        # =========================
        for interaction_type in interactions:

            for interaction in interaction_type:

                try:

                    resnr = interaction.findtext(
                        "resnr"
                    )
                    restype = interaction.findtext(
                        "restype"
                    )
                    reschain = interaction.findtext(
                        "reschain"
                    )

                    resnr_lig = interaction.findtext(
                        "resnr_lig"
                    )
                    restype_lig = interaction.findtext(
                        "restype_lig"
                    )
                    reschain_lig = interaction.findtext(
                        "reschain_lig"
                    )

                    # =====================
                    # METAL COMPLEX
                    # =====================
                    if (
                        interaction_type.tag
                        == "metal_complexes"
                    ):

                        metal_type = (
                            interaction.findtext(
                                "metal_type"
                            )
                        )

                        target_type = (
                            interaction.findtext(
                                "target_type"
                            )
                        )

                        signature = (
                            f"metal_complex:"
                            f"{restype}{resnr}{reschain}"
                            f"->"
                            f"{restype_lig}"
                            f"{resnr_lig}"
                            f"{reschain_lig}"
                            f":metal={metal_type}"
                            f":atom={target_type}"
                        )

                    # =====================
                    # OTHER INTERACTIONS
                    # =====================
                    else:

                        itype = (
                            interaction_type.tag[:-1]
                        )

                        signature = (
                            f"{itype}:"
                            f"{restype}{resnr}{reschain}"
                            f"->"
                            f"{restype_lig}"
                            f"{resnr_lig}"
                            f"{reschain_lig}"
                        )

                    this_frame.add(signature)

                except Exception as e:
                    print(
                        f"Error reading interaction: {e}"
                    )

    frame_interactions[
        frame_id
    ] = this_frame

    # =========================
    # CLEANUP
    # =========================
    os.system(f'rm -f "{pdb_file}"')
    os.system(f'rm -f "{report_path}"')
    os.system("rm -f *_protonated.pdb")
    os.system("rm -f plipfixed*.pdb")

# =========================
# COUNT INTERACTIONS
# =========================
all_interactions = [
    interaction
    for frame in frame_interactions.values()
    for interaction in frame
]

interaction_counts = Counter(
    all_interactions
)

# =========================
# WRITE CSV
# =========================
with open(outfile_name, "w") as f:

    f.write(
        "Interaction_Signature,"
        "Frames_Occurred\n"
    )

    for interaction, count in sorted(
        interaction_counts.items(),
        key=lambda x: x[1],
        reverse=True
    ):

        print(
            f"{interaction}: "
            f"{count} frames"
        )

        f.write(
            f"{interaction},{count}\n"
        )

print(
    f"\nSaved interaction summary to:"
    f" {outfile_name}"
)

# =========================
# TOTAL TIME
# =========================
end_time = time.time()

total_time = end_time - start_time

minutes = int(total_time // 60)
seconds = int(total_time % 60)

print(
    f"\nTotal time taken: "
    f"{minutes} min "
    f"{seconds} sec"
)
