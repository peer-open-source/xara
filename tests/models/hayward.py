from deck_sections import deck_sections
from bent_sections import bent_sections

sections = {
    "deck": deck_sections,
    "bent": bent_sections
}

class Hayward:
    def get_section(self, type, bent, node=None):
        lib = sections[type]
        if bent < 6:
            section = lib["Abut1-Bent6"]

