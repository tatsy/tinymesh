from __future__ import print_function

import os
import sys
import requests
import zipfile
import tarfile
from collections import OrderedDict

class Item:
    def __init__(self, url, name):
        self.url = url
        self.name = name

    def __repr__(self):
        return self.url

target_folder = 'large_meshes'

url_list = OrderedDict([
    ('happy_budda', Item('http://graphics.stanford.edu/pub/3Dscanrep/happy/happy_recon.tar.gz', 'happy_vrip.ply')),
    ('dragon', Item('http://graphics.stanford.edu/pub/3Dscanrep/dragon/dragon_recon.tar.gz', 'dragon_vrip.ply')),
    ('lucy', Item('http://graphics.stanford.edu/data/3Dscanrep/lucy.tar.gz', 'lucy.ply'))
])

def get_file_from_url(url):
    temp = url.split('/')[-1]
    length = temp.rfind('?')
    if length == -1:
        return temp
    else:
        return temp[:length]

def bar_string(rate):
    length = 40
    if rate == 100.0:
        return '=' * length
    else:
        num = int(rate / (100.0 / length))
        return ('=' * num) + '>' + (' ' * (length - num - 1))

def download(url):
    print('Download: {:s}'.format(url))

    filename = get_file_from_url(url)
    r = requests.get(url, stream=True)

    size_total = int(r.headers['content-length'])
    size_get = 0

    chunk_size = 1024
    with open(filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)
                f.flush()
                size_get += chunk_size

                rate = min(100.0, 100.0 * size_get / size_total)
                print('[ %6.2f %% ] [ %s ]' % (rate, bar_string(rate)), end='\r')
        print('\nDownload finished!!')
        return filename

    return None

def unarchive(filename, target, extract_to):
    _, ext = os.path.splitext(filename)
    try:
        if ext == '.gz' or ext == '.tar':
            tar = tarfile.open(filename)
            members = list(filter(lambda mem: target in mem.name, tar.getmembers()))
            if len(members) != 1:
                raise Exception('Multiple files with name "{:s}" are detected!'.format(target))

            with open(extract_to, 'wb') as f:
                f.write(tar.extractfile(members[0]).read())
            tar.close()

        elif ext == '.zip':
            zip_file = zipfile.ZipFile(filename)
            members = list(filter(lambda mem: target in mem.filename, zip_file.infolist()))
            if len(members) != 1:
                raise Exception('Multiple files with name "{:s}" are detected!'.format(target))

            with open(extract_to, 'wb') as f:
                f.write(zip_file.open(members[0].filename).read())
            zip_file.close()

        else:
            raise Exception('Unknown extension: %s' % ext)

    except e as Exception:
        raise e

    print('File is unarchived: {:s}'.format(extract_to))


def process(name, item):
    # Download
    filename = download(item.url)
    # Extract target file
    _, ext = os.path.splitext(item.name)
    extract_to = os.path.join(target_folder, name) + ext
    unarchive(filename, item.name, os.path.join(target_folder, name) + ext)
    # Remove urchive file
    os.remove(filename)

def main():
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    for name, item in url_list.items():
        process(name, item)

if __name__ == '__main__':
    main()