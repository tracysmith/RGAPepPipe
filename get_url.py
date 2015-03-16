#!/usr/bin/python

import sys, argparse, os
from subprocess import call
from multiprocessing.dummy import Pool as ThreadPool

###################################################################
#This is a phython script to download fastq files from ENA
#You can use this directly with the enaFileParser output
###################################################################

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
        os.path.abspath(os.path.expanduser(values)))

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Download fastq files from ENA')
    parser.add_argument("urlFile", help="ERPXXXXXX_download.txt generated from enaFileParser.py", action=FullPaths,
        type=is_file)
    parser.add_argument("-t", "--threads",
        help="Number of threads to use (default: 1)",
        type=int, default=1)
    return parser.parse_args()

def make_urlList(urlFile):
    urls = []
    with open(urlFile, 'r') as infile:
        for line in infile:
            line=line.strip()
            urls.append(line)
    return urls

def download_url(url):
    call('wget {url}'.format(url=url), shell=True)
    ftp = url.split("/")
    index = len(ftp)-1
    filename = ftp[index]
    call('gunzip {filename}'.format(filename=filename), shell=True)

args = get_args()
urls = make_urlList(args.urlFile)

#Make the Pool of workers
pool = ThreadPool(args.threads)
#Open the urls in their own threads and return the results
pool.map(download_url, urls)
pool.close()
pool.join()
