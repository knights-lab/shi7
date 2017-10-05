from urllib.request import urlopen
from sys import platform as _platform
import os


def download_txt_url(path_to_file, url, chunk_size=2 ** 14):
    with urlopen(url) as stream:
        with open(path_to_file, 'wb') as outfile:
            while True:
                chunk = stream.read(chunk_size)
                if not chunk:
                    break
                outfile.write(chunk)


if _platform == "linux" or _platform == "linux2":
    bin_shi7en = 'ninja_shi7_linux'
    bin_split = 'gotta_split_s_linux'
elif _platform == "darwin":
    bin_shi7en = 'ninja_shi7_mac'
    bin_split = 'gotta_split_mac'
elif _platform == "win32":
    bin_shi7en = 'ninja_shi7.exe'
    bin_split = 'gotta_split_s.exe'
else:
    print('Undefined OS\n')
    raise SystemExit()

old_shi7en_bin = 'https://github.com/knights-lab/shi7en/releases/download/0.0.1/' + bin_shi7en
new_shi7en_bin = 'shi7en_trimmer'
old_split_bin = 'https://github.com/knights-lab/shi7en/releases/download/0.0.1/' + bin_split
new_split_bin = 'gotta_split'
download_txt_url(new_shi7en_bin, old_shi7en_bin)
download_txt_url(new_split_bin, old_split_bin)
os.chmod(new_shi7en_bin, 0o755)
os.chmod(new_split_bin, 0o755)

# TODO: maybe add $PATH
