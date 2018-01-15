#!/usr/bin/env python3

"""
Tests for exfi.io
"""


from unittest import TestCase, main

from exfi.io import _coordinate_str_to_tuple

class TestCoordinatesToVariables(TestCase):
    """_coordinate_str_to_tuple(coordinates):
    (string) -> (string, int, int)
    """
    def test_empty_string(self):
        """_coordinate_str_to_tuple: case empty"""
        with self.assertRaises(ValueError):
            _coordinate_str_to_tuple("")

    def test_incorrect_string(self):
        """_coordinate_str_to_tuple: case incorrect"""
        with self.assertRaises(IndexError):
            _coordinate_str_to_tuple("ENSDART00000161035:1")

    def test_correct_string(self):
        """_coordinate_str_to_tuple: case correct"""
        self.assertEqual(
            _coordinate_str_to_tuple("ENSDART00000161035:1-15"),
            ("ENSDART00000161035", 1, 15)
        )

    def test_messy_string(self):
        """_coordinate_str_to_tuple: messy case"""
        self.assertEqual(
            _coordinate_str_to_tuple("TRINITY_g14_c15_i5:1:2-15-12:1-15"),
            ("TRINITY_g14_c15_i5:1:2-15-12", 1, 15)
        )





if __name__ == '__main__':
    main()
