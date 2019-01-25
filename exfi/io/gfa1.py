#!/usr/bin/env python3

"""exfi.io.gfa1.py: submodule for auxliary variables for gfa"""

HEADER_COLS = ['RecordType', 'VersionNumber']
SEGMENT_COLS = ['RecordType', "Name", "Sequence", 'Length']
LINK_COLS = ['RecordType', "From", "FromOrient", "To", "ToOrient", "Overlap"]
CONTAINMENT_COLS = [
    'RecordType', 'Container', 'ContainerOrient', 'Contained',
    'ContainedOrient', 'Pos', 'Overlap'
]
PATH_COLS = ['RecordType', 'PathName', 'SegmentNames', 'Overlaps']
