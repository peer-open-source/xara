from opensees.units.core import UnitManager



systems = ["english_engineering", "si"]


for system in systems:
    units = UnitManager(system)

    assert 1.0 * units.foot == 0.3048 * units.meter
    
    assert 1.0 * units.inch == 0.0245 * units.meter


