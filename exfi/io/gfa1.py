#!/usr/bin/env python3

'''exfi.io.gfa1.py: submodule for auxliary variables for gfa'''

HEADER_COLS = ['record_type', 'version_number']
SEGMENT_COLS = ['record_type', 'name', 'sequence']
LINK_COLS = ['record_type', 'from', 'from_orient', 'to', 'to_orient', 'overlap']
CONTAINMENT_COLS = [
    'record_type', 'container', 'container_orient', 'contained',
    'contained_orient', 'pos', 'overlap'
]
PATH_COLS = ['record_type', 'path_name', 'segment_names', 'overlaps']


HEADER_DTYPES = {'record_type': object, 'version_number': object}
SEGMENT_DTYPES = {
    'record_type': object, 'name': object, 'sequence': object
}
LINK_DTYPES = {
    'record_type': object, 'from': object, 'from_orient': object, 'to': object,
    'to_orient': object, 'overlap': object}
CONTAINMENT_DTYPES = {
    'record_type': object, 'container': object, 'container_orient': object,
    'contained': object, 'contained_orient': object, 'pos': int,
    'overlap': object
}
PATH_DTYPES = {
    'record_type': object, 'path_name': object, 'segment_names': object,
    'overlaps': object
}
