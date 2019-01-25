VERBOSE = "VERBOSE"
DEBUG = "DEBUG"
INFO = "INFO"
WARNING = "WARNING"
ERROR = "ERROR"
FORCED = "FORCED"

# SIDE_CHANNEL is for hidden data that isn't used anywhere else in computation.
SIDE_CHANNEL = "SIDE_CHANNEL"


def log(level, message):
    if level not in [VERBOSE, SIDE_CHANNEL, DEBUG]:
        print(level[0], message)


def assert_eq(a, b, message=""):
    if message: message = message + ": "
    assert a == b, "%sExpected equality of:\n\t%s\n\t%s\n" % (message, a, b)


class FinishNow(Exception):
    pass
