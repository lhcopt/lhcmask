import os
import shutil

folder_list = [
    '../',
    '../../../examples/hl_lhc_collision/',
    '../../../examples/hl_lhc_collision_nobb_b4']

for sixtrack_input_folder in folder_list:
    print(sixtrack_input_folder)
    try:
        for iff in [2,8,16,34]:
            os.system(f"rm {sixtrack_input_folder}/fort.{iff}")
            try:
                shutil.copy(sixtrack_input_folder + f"/fc.{iff}",
                    sixtrack_input_folder + f"/fort.{iff}")
            except Exception:
                print(f"/fc.{iff} not found!")

        with open(sixtrack_input_folder + "/fort.3", "w") as fout:
            with open("fort_parts/fort_beginning.3", "r") as fid_fort3b:
                fout.write(fid_fort3b.read())
            with open(sixtrack_input_folder + "/fc.3", "r") as fid_fc3:
                fout.write(fid_fc3.read())
            with open("fort_parts/fort_end.3", "r") as fid_fort3e:
                fout.write(fid_fort3e.read())
    except Exception as e:
        print('Skipped...')
        print(e)

print('Ended.')
