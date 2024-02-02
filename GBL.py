import os
import sys


linkFormat = "https://genometest2.gs.washington.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position={}%3A{}-{}&hgsid=102808_SZG2N9RA9gcudsoOW98wGyWwauEL"


class Link:
    def __init__(self, pos=None, TAB=None):
        if pos is not None:
            self.input = pos
            self.parsePos()
        elif tab is not None:
            self.input = tab
            self.parseTab()

        self.getLink()

    def parsePos(self):
        s = self.input.split(":")
        self.chr = s[0]
        ss = s[1].split("-")
        self.start = int(ss[0])
        self.end = int(ss[1])

    def parseTab(self):
        s = self.input.split("\t")
        self.chr = s[0]
        self.start = int(s[1])
        self.end = int(s[2])

    def getLink(self):
        self.link = linkFormat.format(self.chr, self.start, self.end)

    def __str__(self):
        return self.link
