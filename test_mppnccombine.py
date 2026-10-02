#!/usr/bin/env python
"""
Copyright 2018 

author:  <scott.wales@unimelb.edu.au>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
from __future__ import print_function
import xarray
import subprocess
import numpy as np
import netCDF4
import pytest
import os
import glob
import sys
import numpy.testing

skip_travis = pytest.mark.skipif(os.environ.get('TRAVIS', None) == 'true', reason="Disabled on Travis")

def run_nccopy(options, infiles, outdir):
    try:
        for f in infiles:
            fout = os.path.join(outdir,os.path.basename(f))
            s = subprocess.check_output(
                    ['nccopy'] + options + [f] + [fout],
                    stderr=subprocess.STDOUT)
            print(s.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(e.stdout.decode('utf-8'))
        raise

# Run the collation program
def run_collate(inputs, output, np=2, args=[]):
    try:
        command = ['mpirun', '--mca', 'btl', 'self,tcp', '--oversubscribe', '-n', '%d'%np, 'mppnccombine-fast', '-o', str(output)] + inputs + args
        print(' '.join(['%s'%x for x in command]))
        subprocess.run(
            command,
            stdout=sys.stdout,
            stderr=subprocess.STDOUT,
            check=True)
    except subprocess.CalledProcessError as e:
        raise
    return xarray.open_dataset(str(output), engine='netcdf4', decode_times=False)

def split_dataset(data, split):
    out = []

    for start in range(0, len(data['x']), split['x']):
        d = data.isel(**{'x': slice(start, start+split['x'])})
        d['x'].attrs['domain_decomposition'] = [1,len(data['x']), 1+start, min(start+split['x'], len(data['x']))]

        out.append(d)

    return out


def split_file(tmpdir, data, split, fill_value=None):
    i = 0
    infiles = []
    for start in range(0, len(data['x']), split['x']):
        infiles.append(str(tmpdir.join('in.%04d.nc'%i)))

        d = data.isel(**{'x': slice(start, start+split['x'])})
        d['x'].attrs['domain_decomposition'] = [1,len(data['x']), 1+start, min(start+split['x'], len(data['x']))]

        d.to_netcdf(infiles[-1],
                encoding = {'a': {
                    'chunksizes': (2,),
                    'zlib': True,
                    'shuffle': True,
                    'complevel': 4,
                    '_FillValue': fill_value,
                    },
                    }
                )
        d.close()
        i+=1

    return infiles


def test_simple(tmpdir):
    inpath = str(tmpdir.join('test.nc'))
    outpath = tmpdir.join('out.nc')
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(5))
            },
            coords = {
                'x': np.arange(5),
            })

    d.to_netcdf(inpath, encoding={
            'a': {
                'chunksizes': (2,),
                'zlib': True,
                'shuffle': True,
                'complevel': 4,
                },
        })

    c = run_collate([inpath], outpath)

    assert (c.a.data == d.a.data).all()
    numpy.testing.assert_array_equal(c.a.data, d.a.data)


def test_split_on_boundary(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4))
            },
            coords = {
                'x': np.arange(4),
            })

    infiles = split_file(tmpdir, d, {'x': 2})

    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath)

    assert (c.a.data == d.a.data).all()
    assert (c.x.data == d.x.data).all()
    
    assert 'domain_decomposition' not in c.x.attrs


def test_split_off_boundary(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(6))
            },
            coords = {
                'x': np.arange(6),
            })

    infiles = split_file(tmpdir, d, {'x': 3})

    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath)

    assert (c.a.data == d.a.data).all()
    assert (c.x.data == d.x.data).all()
    
@pytest.mark.parametrize('fill_value, expected', [
    (None, 0.0),
    (-1.0e20, -1.0e20),
])
def test_missing_tile(tmpdir, fill_value, expected):
    """
    A tile eliminated by a mask table is never written, so its region of the
    collated output is covered by no input file. Those cells should take the
    declared _FillValue, or zero if none is declared.
    """
    d = xarray.Dataset(
            {
                'a': (['x'], np.ones(6))
            },
            coords = {
                'x': np.arange(6),
            })

    infiles = split_file(tmpdir, d, {'x': 2}, fill_value=fill_value)

    # Drop the middle tile
    del infiles[1]

    outpath = tmpdir.join('out.nc')
    run_collate(infiles, outpath)

    with netCDF4.Dataset(str(outpath)) as nc:
        a = nc.variables['a']
        a.set_auto_mask(False)
        numpy.testing.assert_array_equal(a[:], [1, 1, expected, expected, 1, 1])

        if fill_value is None:
            assert '_FillValue' not in a.ncattrs()
        else:
            assert a.getncattr('_FillValue') == fill_value


def test_different_compression(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4))
            },
            coords = {
                'x': np.arange(4),
            })

    split = split_dataset(d, {'x': 2})

    infiles = [str(x) for x in [tmpdir.join('000.nc'), tmpdir.join('001.nc')]]

    # Save with different compression settings
    split[0].to_netcdf(infiles[0],
            encoding = {'a': {
                'chunksizes': (2,),
                'zlib': True,
                'shuffle': True,
                'complevel': 4,
                }})
    split[1].to_netcdf(infiles[1],
            encoding = {'a': {
                'chunksizes': (2,),
                'zlib': True,
                'shuffle': False,
                'complevel': 6,
                }})

    outpath = tmpdir.join('out.nc')

    c = run_collate(infiles, outpath)

    np.testing.assert_array_equal(d.a, c.a)
    np.testing.assert_array_equal(d.x, c.x)


def test_compression_override(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4))
            },
            coords = {
                'x': np.arange(4),
            })

    infiles = split_file(tmpdir, d, {'x': 2})

    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath, args=['-d','8','--no-shuffle'])

    assert c.a.encoding['complevel'] == 8
    assert c.a.encoding['shuffle'] == False
    np.testing.assert_array_equal(d.a, c.a)

def test_clobber(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4))
            },
            coords = {
                'x': np.arange(4),
            })

    infiles = split_file(tmpdir, d, {'x': 2})

    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath)
    c.close()

    with pytest.raises(subprocess.CalledProcessError):
        c = run_collate(infiles, outpath)

    c = run_collate(infiles, outpath, args=['-f'])

def test_clean(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4))
            },
            coords = {
                'x': np.arange(4),
            })

    infiles = split_file(tmpdir, d, {'x': 2})

    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath, args=['-r'])

    assert not os.path.isfile(infiles[0])
    assert not os.path.isfile(infiles[1])
    assert os.path.isfile(outpath)

def test_1degree(tmpdir):

    # This test dataset is masked, has a missing tile, as there is
    # no data in that tile, and has inconsistent chunk sizes
    testdir = 'test_data'
    testdatadir = os.path.join(testdir, 'output_1deg_masked')

    infiles = glob.glob(os.path.join(testdatadir,'ocean_month.nc.????'))

    # Open file collated by mppnccombine
    d = xarray.open_dataset(os.path.join(testdatadir,'ocean_month.nc'),decode_times=False)

    # No compression
    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath)

    assert d.equals(c)

    # Compression + shuffle
    outpath = tmpdir.join('out_shuff.nc')
    c = run_collate(infiles, outpath, args=['--shuffle','--deflate','5'])

    assert d.equals(c)

def test_1degree_nc3(tmpdir):

    # This test dataset is masked, has a missing tile, as there is
    # no data in that tile, and has inconsistent chunk sizes
    testdir = 'test_data'
    testdatadir = os.path.join(testdir, 'output_1deg_masked')

    infiles = glob.glob(os.path.join(testdatadir,'ocean_month.nc.????'))

    # Open file collated by mppnccombine
    d = xarray.open_dataset(os.path.join(testdatadir,'ocean_month.nc'),decode_times=False)

    # Convert files to netCDF Classic format and test
    testdir = os.path.join(tmpdir,testdatadir+"_nc3")
    os.makedirs(testdir)
    run_nccopy(['-3'],infiles,testdir)

    infiles = glob.glob(os.path.join(testdir,'ocean_month.nc.????'))

    # No compression
    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath)

    assert d.equals(c)

    # Compression + shuffle
    outpath = tmpdir.join('out_shuff.nc')
    c = run_collate(infiles, outpath, args=['--shuffle','--deflate','5'])

    assert d.equals(c)
    assert c.sst.encoding['complevel'] == 5
    assert c.sst.encoding['shuffle'] == True

def test_unlimited(tmpdir):
    inpath = str(tmpdir.join('test.nc'))
    outpath = tmpdir.join('out1.nc')
    d = xarray.Dataset(
            {
                'a': (['time','x'], np.random.rand(6,5))
            },
            coords = {
                'x': np.arange(5),
                'time': np.arange(6),
            })

    d.to_netcdf(inpath,
            unlimited_dims = ['time'],
            encoding={
            'a': {
                'chunksizes': (2,2),
                'zlib': True,
                'shuffle': True,
                'complevel': 4,
                },
        })
    c = run_collate([inpath], outpath)
    assert c.encoding['unlimited_dims'] == set(['time'])

    outpath = tmpdir.join('out2.nc')
    d = xarray.Dataset(
            {
                'a': (['time','x'], np.random.rand(1,5))
            },
            coords = {
                'x': np.arange(5),
                'time': np.arange(1),
            })

    d.to_netcdf(inpath,
            unlimited_dims = ['time'],
            encoding={
            'a': {
                'chunksizes': (1,2),
                'zlib': True,
                'shuffle': True,
                'complevel': 4,
                },
        })
    c = run_collate([inpath], outpath, args=['--force'])
    assert c.encoding['unlimited_dims'] == set(['time'])

def test_many_files(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4000))
            },
            coords = {
                'x': np.arange(4000),
            })

    infiles = split_file(tmpdir, d, {'x': 2})

    outpath = tmpdir.join('out.nc')
    c = run_collate([tmpdir.join('*.nc')], outpath)

    assert d.equals(c)

def test_min_dim(tmpdir):
    """
    Minimum domain value should be respected 
    """
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4))
            },
            coords = {
                'x': np.arange(4),
            })
    d['x'].attrs['domain_decomposition'] = [101, 104, 101, 104]

    infile = str(tmpdir.join('offset.nc'))
    d.to_netcdf(infile, encoding={'a':{'chunksizes':(2,)}})

    outpath = tmpdir.join('out.nc')
    c = run_collate([infile], outpath)

    assert c.x.size == 4

def test_partial_vertical_chunk(tmpdir):
    inpath = str(tmpdir.join('test.nc'))
    outpath = tmpdir.join('out.nc')
    d = xarray.Dataset(
            {
                'a': (['z','x'], numpy.zeros((3,1)))
            },
            coords = {
                'x': np.arange(1),
                'z': np.arange(3),
            })
    d.x.attrs['domain_decomposition'] = [1, 1, 1, 1]

    d.to_netcdf(inpath, encoding={
            'a': {
                'chunksizes': (2,1),
                },
        })

    c = run_collate([inpath], outpath)

    numpy.testing.assert_array_equal(c.a.data, d.a.data)

@skip_travis
def test_four_processes(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(8))
            },
            coords = {
                'x': np.arange(8),
            })

    infiles = split_file(tmpdir, d, {'x': 2})

    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath, np=4)

    numpy.testing.assert_array_equal(c.a.data, d.a.data)

def test_no_match(tmpdir):
    """
    Error if the input file is missing
    """
    infile = str(tmpdir.join('missing.nc'))
    outpath = tmpdir.join('out.nc')

    with pytest.raises(subprocess.CalledProcessError):
        c = run_collate([infile], outpath)
