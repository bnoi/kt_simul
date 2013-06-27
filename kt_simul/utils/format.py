import collections

def pretty_dict(indict, indent_level=0):
    """
    Return a pretty to display string from an ordered dict.
    "separate" key put a new line in the output.
    """

    out = ""
    indentation = " " * 4 * indent_level
    for key, value in indict.items():
        if key == "separate":
            out += "\n"
        elif isinstance(value, collections.OrderedDict):
            value_out = pretty_dict(value, indent_level + 1)
            out += "%s%s : \n%s\n" % (indentation, key, value_out)
        else:
            out += "%s%s : %s\n" % (indentation, key, value)

    return out
