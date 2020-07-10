import warnings


class DesignError(Exception):
    pass


def deprecation(message):
    warnings.warn(message, stacklevel=3)


class EquationWarning(Warning):
    pass