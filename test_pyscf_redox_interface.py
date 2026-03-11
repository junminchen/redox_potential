from pyscf_redox_interface import (
    ChargeStateSpec,
    DFTStateResult,
    build_redox_config_entry,
    load_geometry,
)


def test_load_pdb_geometry():
    atoms = load_geometry("pdb_bank/EC.pdb", "pdb")
    assert len(atoms) > 0
    assert atoms[0][0].isalpha()


def test_build_redox_config_entry():
    states = [
        ChargeStateSpec(state=0, charge=0, spin=0),
        ChargeStateSpec(state=-1, charge=-1, spin=1),
    ]
    results = [
        DFTStateResult(
            state=0,
            charge=0,
            spin=0,
            total_energy_hartree=-100.0,
            relative_energy_kjmol=0.0,
            effective_offset_kjmol=0.0,
            predicted_half_wave_v_vs_vacuum=None,
            predicted_half_wave_v_vs_li=None,
        ),
        DFTStateResult(
            state=-1,
            charge=-1,
            spin=1,
            total_energy_hartree=-100.1,
            relative_energy_kjmol=-10.0,
            effective_offset_kjmol=25.0,
            predicted_half_wave_v_vs_vacuum=-0.259,
            predicted_half_wave_v_vs_li=1.141,
        ),
    ]
    entry = build_redox_config_entry(
        molecule_name="EC",
        residue_names=["EC"],
        states=states,
        state_results=results,
    )
    assert entry["state_free_energy_offsets_kjmol"]["0"] == 0.0
    assert entry["state_free_energy_offsets_kjmol"]["-1"] == 25.0
    assert "0,-1" in entry["electron_affinities"]
