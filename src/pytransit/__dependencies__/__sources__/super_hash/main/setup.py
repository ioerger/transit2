import setuptools
import toml
import os

# 
# get the data out of the toml file
# 
toml_info = toml.load("../pyproject.toml")
package_info = {**toml_info["tool"]["poetry"], **toml_info["tool"]["extra"]}

# 
# file system helpers
# 
def final_target_of(path):
    # resolve symlinks
    if os.path.islink(path):
        have_seen = set()
        while os.path.islink(path):
            path = os.readlink(path)
            if path in have_seen:
                return None # circular broken link
            have_seen.add(path)
    return path

def is_folder(path):
    # resolve symlinks
    final_target = final_target_of(path)
    if not final_target:
        return False
    return os.path.isdir(final_target)

def iterate_folder_paths_in(path, recursively=False, have_seen=None):
    if recursively and have_seen is None: have_seen = set()
    if is_folder(path):
        with os.scandir(path) as iterator:
            for entry in iterator:
                entry_is_folder = False
                if entry.is_dir():
                    entry_is_folder = True
                    each_subpath_final_target = os.path.join(path, entry.name)
                # check if symlink to dir
                elif entry.is_symlink():
                    each_subpath = os.path.join(path, entry.name)
                    each_subpath_final_target = final_target_of(each_subpath)
                    if is_folder(each_subpath_final_target):
                        entry_is_folder = True
                
                if entry_is_folder:
                    yield os.path.join(path, entry.name)
                    if recursively:
                        have_not_already_seen_this = each_subpath_final_target not in have_seen
                        have_seen.add(each_subpath_final_target) # need to add it before recursing
                        if have_not_already_seen_this:
                            yield from iterate_folder_paths_in(each_subpath_final_target, recursively, have_seen)


# 
# get the data out of the readme file
# 
with open("../README.md", "r") as file_handle:
    long_description = file_handle.read()

# 
# generate the project
#  
setuptools.setup(
    name=package_info["name"],
    version=package_info["version"],
    description=package_info["description"],
    url=package_info["url"],
    author=package_info["author"],
    author_email=package_info["author_email"],
    license=package_info["license"],
    packages=[package_info["name"]],
    install_requires=[
        # examples:
        # 'aiohttp >= 3.7.4',
        # 'requests == 2.26.0',
        # <, >, <=, >=, == or !=
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Programming Language :: Python',
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    long_description=long_description,
    long_description_content_type="text/markdown",
    python_requires='>=3.6',
    package_data={
        package_info["name"]: [
            (each+"/*.py")[ len(package_info["name"])+1 : ] 
                for each in iterate_folder_paths_in(package_info["name"], recursively=True)
        ],
    },
)