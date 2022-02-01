#!/usr/bin/env python3

class Donor:
    def __init__(self, donor_id, age, sex, cells=[]):
        self.donor_id = donor_id
        self.age = age
        self.sex = sex
        self.cells = cells

class Cell:
    def __init__(self, cell_barcode, sequenced_areas = []):
        self.cell_barcode = cell_barcode
        self.sequenced_areas = sequenced_areas

class SequencedArea:
    def __init__(self, umis = []):
        self.umis = umis

class UMI:
    def __init__(self, UMI_id, read_sets = []):
        self.UMI_id = UMI_id
        self.read_sets = read_sets

class ReadSet:
    def __init__(self, reads = []):
        self.reads = reads

class Read:
    def __init__(self, sequence):
        self.sequence = sequence

