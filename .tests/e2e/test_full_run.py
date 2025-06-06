import os
import pytest
import shutil
import subprocess as sp
from pathlib import Path


@pytest.fixture
def setup(tmp_path):
    reads_fp = Path(".tests/data/reads/").resolve()
    ref_fp = Path(".tests/data/ref/").resolve()
    mapping_fp = Path(".tests/data/mapping.csv").resolve()
    blast_db_fp = tmp_path / "blast_db"
    blast_db_fp.mkdir(parents=True, exist_ok=True)
    (blast_db_fp / "db.faa").touch()
    phold_db_fp = tmp_path / "phold_db"
    phold_db_fp.mkdir(parents=True, exist_ok=True)

    project_dir = tmp_path / "project/"

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = project_dir / "sunbeam_config.yml"

    config_str = f"sbx_induceome: {{reference_fp: {ref_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"sbx_induceome: {{mapping_fp: {mapping_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"sbx_induceome: {{blast_db: {blast_db_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"sbx_induceome: {{phold_db: {phold_db_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield tmp_path, project_dir

    shutil.rmtree(tmp_path)


@pytest.fixture
def run_sunbeam(setup):
    tmp_path, project_dir = setup
    output_fp = project_dir / "sunbeam_output"
    log_fp = output_fp / "logs"
    stats_fp = project_dir / "stats"

    sbx_proc = sp.run(
        [
            "sunbeam",
            "run",
            "--profile",
            project_dir,
            "all_induceome",
            "--directory",
            tmp_path,
            "-n",
        ],
        capture_output=True,
        text=True,
    )

    print("STDOUT: ", sbx_proc.stdout)
    print("STDERR: ", sbx_proc.stderr)

    if os.getenv("GITHUB_ACTIONS") == "true":
        try:
            shutil.copytree(log_fp, "logs/")
            shutil.copytree(stats_fp, "stats/")
        except FileNotFoundError:
            print("No logs or stats directory found.")

    output_fp = project_dir / "sunbeam_output"
    benchmarks_fp = project_dir / "stats/"

    yield output_fp, benchmarks_fp, sbx_proc


def test_full_run(run_sunbeam):
    output_fp, benchmarks_fp, proc = run_sunbeam

    assert proc.returncode == 0, f"Sunbeam run failed with error: {proc.stderr}"
