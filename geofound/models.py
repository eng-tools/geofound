from sfsimodels import models


class Soil(models.Soil):
    gwl = 0.0  # ground water level


class Foundation(models.Foundation):
    q_ult = 0.0

    def __repr__(self):
        return self.length, self.width, self.height
