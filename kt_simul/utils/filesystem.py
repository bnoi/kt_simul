"""
"""

import os
import glob

def tree(path):
    """
    """

    all_files = []
    for root, subdir, files in os.walk(path):
        for f in files:
            relpath = os.path.relpath(os.path.join(root, f), path)
            all_files.append(relpath)

    return all_files