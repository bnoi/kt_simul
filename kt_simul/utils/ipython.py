def in_ipython():
    try:
        __IPYTHON__
    except NameError:
        return False
    else:
        return True
