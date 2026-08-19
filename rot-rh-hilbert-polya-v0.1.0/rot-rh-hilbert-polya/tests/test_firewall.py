from pathlib import Path

def test_v100_imports_mpmath_only_after_hard_freeze_section():
    root = Path(__file__).resolve().parents[1]
    text = (root/'research/legacy_pipeline/rot_rh_full_arithmetic_tail_v100_1.py').read_text(encoding='utf-8')
    hard = text.index('# HARD FREEZE')
    saved = text.index('np.savez_compressed', hard)
    target = text.index('import mpmath as mp')
    assert hard < saved < target


def test_v100_declares_no_target_data_pre_freeze():
    root = Path(__file__).resolve().parents[1]
    text = (root/'research/legacy_pipeline/rot_rh_full_arithmetic_tail_v100_1.py').read_text(encoding='utf-8')
    assert 'Xi/zeta/known zeros pre-freeze: NONE' in text
    assert '"Xi_used_pre_freeze":False' in text
    assert '"known_zeros_used_pre_freeze":False' in text
