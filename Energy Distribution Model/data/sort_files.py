import os
import shutil

def organize_files_by_filename(source_dir, target_dirs_and_strings):
    """
    Organizes files in `source_dir` into folders based on substrings in their filenames.

    Parameters:
        source_dir (str): Path to the directory containing files to organize.
        target_dirs_and_strings (list of tuples): Each tuple contains:
            - target_dir (str): Path to the target directory where matching files will be moved.
            - search_string (str): Substring to search for in the filenames.
    """
    # Ensure the source directory exists
    if not os.path.exists(source_dir):
        print(f"Source directory '{source_dir}' does not exist.")
        return

    # Create target directories if they don't exist
    for target_dir, _ in target_dirs_and_strings:
        os.makedirs(target_dir, exist_ok=True)

    # Iterate over files in the source directory
    for filename in os.listdir(source_dir):
        file_path = os.path.join(source_dir, filename)

        # Skip directories
        if not os.path.isfile(file_path):
            continue

        # Check filenames for each substring
        for target_dir, search_string in target_dirs_and_strings:
            if search_string in filename:
                # Move the file to the corresponding folder
                shutil.move(file_path, os.path.join(target_dir, filename))
                print(f"Moved '{filename}' to '{target_dir}'")
                break  # Stop checking other strings once a match is found

# Example usage
source_directory = "0.95"
target_directories_and_strings = [
    (source_directory + "//lab-energies", "Ek.asc"),
    (source_directory + "//e-densities", "ne.asc"),
    (source_directory + "//data", ".dat"),
]

organize_files_by_filename(source_directory, target_directories_and_strings)
