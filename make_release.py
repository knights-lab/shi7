#!/usr/bin/env python
import tempfile
from distutils.dir_util import copy_tree
import os
import shutil

entry_points = ("shi7", "shi7_learning")

def create_main_template(entry_point):
    template = "#!/usr/bin/env python\nimport shi7\nshi7.{entry_point}.main()\n"
    for entry_point in entry_points:
        yield entry_point, template.format(entry_point=entry_point)

def line_prepender(filename, line):
    with open(filename, 'rb+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip(b'\r\n') + b'\n' + content + b'\n')

for entry_point, main_file in create_main_template(entry_points):
    with tempfile.TemporaryDirectory() as tmpdirname:
        os.makedirs(os.path.join(tmpdirname, 'shi7'))
        copy_tree("shi7", os.path.join(tmpdirname, 'shi7'))
        with open(os.path.join(tmpdirname, '__main__.py'), 'w') as inf:
            inf.write(main_file)
        shutil.make_archive(entry_point, 'zip', tmpdirname)
        # shutil.move(entry_point + '.zip', entry_point + '.pyz')
        # shutil.move(os.path.join(tmpdirname, entry_point),"./")
        # line_prepender(entry_point + '.pyz', b"#!/usr/bin/env python")



