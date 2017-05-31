#!/usr/bin/env python3


def filter_by_length(iterables, length):
    """Filter each of the elements in iterables by length"""
    return (element for element in iterables if len(element.seq) >= length)
