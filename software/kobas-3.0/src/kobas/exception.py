from sys import stderr

class FastaIOError(IOError):
    pass

class ResNulError(Exception):
    pass

class TabformatErr(Exception):
    pass

class FastaErr(Exception):
    pass

class AnnotError(Exception):
    pass

class StatError(ArithmeticError):
    pass

### new for gsea
class GmtformatErr(Exception):
    pass

class GmtdatabaseErr(Exception):
    pass

class GmtvalueErr(Exception):
    pass

class ClsformatErr(Exception):
    pass
class ClsvalueErr(Exception):
    pass

class GctformatErr(Exception):
    pass
class GctvalueErr(Exception):
    pass
class FormatErr(Exception):
    pass

def error(msg):
    err_msg  = '%s: %s\n' % (msg.__class__, msg)
    stderr.write(err_msg)
