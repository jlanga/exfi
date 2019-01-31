#!/usr/bin/env python3

"""exfi.logger.py: submodule to set up the logger"""

import logging


def set_up_logger(args):
    """Get from the parserd argparse the flags related to logging and set up
    the logger
    """

    # Set up logger
    logger = logging.getLogger()
    logging.basicConfig(
        format='%(asctime)s\t%(module)s\t%(message)s',
        level=logging.ERROR
    )
    if args["verbose"]:
        logger.setLevel(logging.INFO)
    if args["debug"]:
        logger.setLevel(logging.DEBUG)

    return logger
