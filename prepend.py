"""
Script to prepend a text to all Python files in the current directory and subdirectories. Used for development purposes.
"""

# Standard Library Imports

# Local Library Imports

####################################################################################################

import os
import glob

# Define the text to prepend
text_to_prepend = '''"""

"""

# Standard Library Imports

# Local Library Imports

####################################################################################################
'''

# Get the current working directory
cwd = os.getcwd()

# Find all Python files in the current directory and subdirectories
for filename in glob.glob(os.path.join(cwd, '**/*.py'), recursive=True):
    with open(filename, 'r+') as f:
        content = f.read()
        if not content.strip():  # Check if the file is empty
            f.seek(0, 0)
            f.write(text_to_prepend + '\n' + content)
