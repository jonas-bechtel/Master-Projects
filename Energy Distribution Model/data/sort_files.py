import os
import shutil

def organize_files_by_filename(source_dir, target_dirs_and_strings, extra_filter=None):
    """
    Organizes files based on substrings in filenames and optional extra filter.

    Parameters:
        source_dir (str): Source directory path.
        target_dirs_and_strings (list of tuples): [(target_dir, search_string), ...]
        extra_filter (str): Optional additional string that must also be in the filename.
    """
    if not os.path.exists(source_dir):
        print(f"Source directory '{source_dir}' does not exist.")
        return

    # Create target directories
    for target_dir, _ in target_dirs_and_strings:
        os.makedirs(target_dir, exist_ok=True)

    for filename in os.listdir(source_dir):
        file_path = os.path.join(source_dir, filename)

        if not os.path.isfile(file_path):
            continue

        for target_dir, search_string in target_dirs_and_strings:
            # Check both conditions: search_string AND extra_filter (if provided)
            if search_string in filename and (extra_filter is None or extra_filter in filename):
                shutil.move(file_path, os.path.join(target_dir, filename))
                print(f"Moved '{filename}' to '{target_dir}'")
                break

# Example usage:
source_directory = "C60\\mbrc_merged_11.3_v2"
target_directories_and_strings = [
    (os.path.join(source_directory, "lab-energies"), "Ek.asc"),
    (os.path.join(source_directory, "e-densities"), "ne.asc"),
    (os.path.join(source_directory, "data"), ".dat"),
]

# Only move files that ALSO contain "run1" in the filename
organize_files_by_filename(source_directory, target_directories_and_strings, extra_filter="")
