This directory must contain templates which specify what to write as common global parameters (at the end of a .model file).
This is an alternative solution to using presets defined by functions inside io_star_common.cpp
To use a template instead of the default, this directory must contain a file named use_template.cfg with two lines:
First line: A boolean that specify if we use the template
Second line: The name of the template to be used. The file is expected into this directory
