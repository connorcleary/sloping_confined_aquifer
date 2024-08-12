import os

def main():
    old_name = "saltelli"
    new_name = "test_saltelli2"

    files = os.listdir(f"/home/superuser/sloping_confined_aquifer/model_files_backup/")

    for file in files:
        if file.startswith(old_name):
            try:
                os.rename(f"/home/superuser/sloping_confined_aquifer/model_files_backup/"+file, f"/home/superuser/sloping_confined_aquifer//model_files_backup/"+file.replace(old_name, new_name))
            except:
                os.rmdir(f"/home/superuser/sloping_confined_aquifer/model_files_backup/"+file)
                os.rename(f"/home/superuser/sloping_confined_aquifer/model_files_backup/"+file, f"/home/superuser/sloping_confined_aquifer//model_files_backup/"+file.replace(old_name, new_name))


if __name__ == "__main__":
    main()