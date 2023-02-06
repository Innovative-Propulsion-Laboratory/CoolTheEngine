# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 11:05:04 2020

@author: julie
"""


class ProgressBar:

    def __init__(self, valmax, maxbar, title):
        self.valmax = valmax  # Maximum value reached by the bar
        self.maxbar = maxbar  # Maximum number of █ in the bar (size of the bar)
        self.title = title  # Title of the bar

    def update(self, val):
        import sys
        # format
        if val > self.valmax:
            val = self.valmax

        # process
        perc = round((float(val) / float(self.valmax)) * 100)  # Percentage computation
        scale = 100.0 / float(self.maxbar)
        bar = int(perc / scale)  # Number of bar needed at each update

        # render 
        out = '\r█ %20s |%s%s| %3d %%  █' % (self.title, '█' * bar, ' ' * (self.maxbar - bar), perc)
        sys.stdout.write(out)

"""
For the out :
Writing of the title on 20 characters (%20s = title) then | 
then filling of the bar (%s %s = '█' * bar, ' ' * (self.maxbar - bar)) 
then | then the pourcentage in three digits (%3d) then the symbol % (%%)
"""