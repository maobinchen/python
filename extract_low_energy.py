#! /usr/bin/env python2.7
# Extracts the low energy k models from a user-specified silent file
# in O(n) time and O(k) space. The resulting models can be extracted
# as PDBs into a user-specified directory or transferred directly to
# another silent file. Specific energy terms can be subtracted from
# the score.
#
# @author Christopher Miles (cmiles@uw.edu)
from heapq import *
import argparse, os, os.path, re, shutil

def extract_silent(options, tags):
    '''Writes the low-scoring entries from the input silent file to the output silent file'''
    hash = set(tags)

    with open(options['silent_out'], 'w') as fout:
        with open(options['silent_in']) as fin:
            fout.write(fin.readline())  # SEQUENCE
            fout.write(fin.readline())  # SCORE header

            for line in fin:
                cols = line.split()
                tag = cols[-1].strip()

                if tag in hash:
                    fout.write(line)


def extract_dir(options, tags):
    '''Writes the low-scoring entries from the input silent file to PDBs in the output directory'''
    bin = options['src'] + '/bin/extract_pdbs.default.' + options['platform'] + options['compiler'] + options['optimization']
    assert os.path.exists(bin), 'Failed to locate extract_pdbs binary'

    extract = [bin,
               '-silent_read_through_errors',
               '-database ' + options['db'],
               '-in:file:silent_struct_type binary',
               '-in:file:silent',
               options['silent_in'],
               '-tags',
               ' '.join(tags)]

    # Execute the command, then move extracted files to the output directory
    os.system(' '.join(extract))

    for tag in tags:
        src = tag + '.pdb'
        dst = options['output_dir'] + '/' + src

        if os.path.exists(src):
            shutil.move(src, dst)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # Input options
    parser.add_argument('--silent_in', '-i', required = True, help = 'Input silent file')

    # Output options
    parser.add_argument('--silent_out', '-o', help = 'Output silent file')
    parser.add_argument('--output_dir', '-d', help = 'Directory to store extracted models')
    parser.add_argument('--num_models', '-n', type = int, default = 100, help = 'Number of models to extract')

    # Filtering options
    parser.add_argument('--exclude_terms', default = '', nargs = '+', help = 'Energy terms to subtract from score. Supports regular expressions.')

    # Rosetta-specific parameters
    parser.add_argument('--src', default = '/work/brunette/src/Rosetta/main/source', help = 'Path to rosetta_source')
    parser.add_argument('--db', default = '/work/brunette/src/Rosetta/main/database', help = 'Path to rosetta_database')
    parser.add_argument('--optimization', default = 'release', help = 'Optimization level')
    parser.add_argument('--platform', default = 'linux', help = 'Platform')
    parser.add_argument('--compiler', default = 'gcc', help = 'Compiler')
    options = vars(parser.parse_args())

    assert (options['silent_out'] or options['output_dir']), 'Must provide either --silent_out or --output_dir'
    assert os.path.exists(options['silent_in']), 'Failed to open input silent file'
    assert options['num_models'] >= 0


    # Removes the output directory if it already exists
    if options['output_dir']:
        if os.path.exists(options['output_dir']):
            shutil.rmtree(options['output_dir'])
        os.mkdir(options['output_dir'])


    with open(options['silent_in']) as file:
        file.readline()

        # Columns containing score, tag
        header = file.readline().split()
        idx_score = header.index('score')
        idx_descr = header.index('description')

        # Columns containing score terms to exclude
        excluded_idxs = []
        for e in options['exclude_terms']:
            pattern = re.compile(e)

            for (i, term) in enumerate(header):
                if pattern.search(term):
                    print 'Excluding %s, index %d' % (term, i)
                    excluded_idxs.append(i)

        heap = []

        for line in file:
            line = line.strip()
            if not line.startswith('SCORE'):
                continue

            cols = line.split()

            try:
                score = float(cols[idx_score])
                descr = cols[idx_descr]

                # Subtract all excluded terms from the score
                additional = [float(cols[i]) for i in excluded_idxs]
                score -= sum(additional)

                # Negate score because we're inserting into a min heap
                entry = (-score, descr)

                if len(heap) < options['num_models']:
                    heappush(heap, entry)
                else:
                    heappushpop(heap, entry)
            except:
                pass

    tags = []
    while heap:
        score, descr = heappop(heap)
        tags.append(descr)

    # Write the low-energy decoys to the appropriate sink
    output_func = extract_silent if options['silent_out'] else extract_dir
    output_func(options, tags)
