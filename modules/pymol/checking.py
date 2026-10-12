cmd = __import__("sys").modules["pymol.cmd"]

def _raising(code=-1, _self=cmd):
    # WARNING: internal routine, subject to change
    return is_error(code)


def is_string(obj):
    return isinstance(obj, (str, bytes))


def is_list(obj):
    return isinstance(obj, list)


def is_dict(obj):
    return isinstance(obj, dict)


def is_tuple(obj):
    return isinstance(obj, tuple)


def is_sequence(obj):
    return isinstance(obj, (list, tuple))


def is_error(result):  # errors are always negative numbers
    if isinstance(result, int):
        return (result < 0)
    return 0


def is_ok(result):  # something other than a negative number
    return not is_error(result)
