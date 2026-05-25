import logging
import platform
import subprocess
import sys
from pathlib import Path

import clapy

from _commons import console, ensure_dir
from program_config import DEFAULT_PROGRAMS_PATHS_FILE

INSTALL_PROGRAMS_COMMAND_NAME = "install-programs"

# fmt: off
cli: clapy.Command = (
    clapy.Command(INSTALL_PROGRAMS_COMMAND_NAME)
        .help("Installs the necessary programs for you.")
        .arguments([
            clapy.Arg("programs-dir")
                .help("Directory to install the programs in.")
                .long("--dir")
                .short("-d")
                .value_parser(Path),
            clapy.Arg("programs-paths-file")
                .help("File to store the installed programs paths, usable through the '--programs-paths' argument.")
                .long("--paths-file")
                .default(DEFAULT_PROGRAMS_PATHS_FILE)
        ])
)
# fmt: on


def run_micromamba(
    logger: logging.Logger,
    target_programs_dir: Path,
    target_programs_paths_file: Path,
) -> None:
    """
    Will install any requried dependencies and configure the program to use them.
    """
    micromamba_bin: str | None = None
    match (platform.system(), platform.architecture()[0]):
        case ("Linux", "64bit"):
            micromamba_bin = "micromamba-linux-64"
        case ("Darwin", "64bit"):
            micromamba_bin = "micromamba-osx-64"

    if micromamba_bin is None:
        logger.critical(f"OS not supported: {platform.system()}")
        sys.exit(1)

    micromamba_path = Path(__file__).parent / "micromamba" / micromamba_bin
    if not micromamba_path.exists():
        logger.critical(
            "Couldn't find the necessary micromamba executable for this OS. Make sure to download the correct release for your OS (https://github.com/DarthPapalo/phagefind/releases)"
        )

    env_lock_path = Path(__file__).parent / "programs_env-lock.yml"

    target_programs_dir_override: Path = target_programs_dir.resolve()

    ensure_dir(
        logger,
        target_programs_dir_override,
        True,
        "can't install programs",
    )

    with console.status("Installing programs..."):
        failed: bool = False
        try:
            subprocess.run(
                [
                    str(micromamba_path),
                    "create",
                    "-p",
                    str(target_programs_dir_override),
                    "-f",
                    str(env_lock_path),
                    "--yes",
                    "--no-rc",
                    "--no-env",
                    "--retry-clean-cache",
                ],
                stderr=subprocess.STDOUT,
                stdout=subprocess.DEVNULL,
                check=True,
            )
        except subprocess.CalledProcessError:
            logger.critical(
                "Couldn't install required programs with micromamba, exiting"
            )
            failed = True
        except PermissionError:
            logger.critical(
                f"Micromamba executable doesn't have execution permission. Try running `chmod +x {micromamba_path}`"
            )
            failed = True
        finally:
            if failed:
                target_programs_dir_override.rmdir()
                sys.exit(1)

    logger.info(f"Programs installed in '{target_programs_dir_override}'")

    with target_programs_paths_file.open("w") as f:
        f.writelines(
            [
                f'samtools = "{str(target_programs_dir_override / "bin" / "samtools")}"\n\n',
                f'fastp = "{str(target_programs_dir_override / "bin" / "fastp")}"\n\n',
                f'spades = "{str(target_programs_dir_override / "bin" / "spades.py")}"\n\n',
                "[blastplus]\n",
                f'blastn = "{str(target_programs_dir_override / "bin" / "blastn")}"\n',
                f'makeblastdb = "{str(target_programs_dir_override / "bin" / "makeblastdb")}"\n\n',
                "[mummer]\n",
                f'nucmer = "{str(target_programs_dir_override / "bin" / "nucmer")}"\n',
                f'show_aligns = "{str(target_programs_dir_override / "bin" / "show-aligns")}"\n',
                f'dnadiff = "{str(target_programs_dir_override / "bin" / "dnadiff")}"\n',
            ]
        )

    logger.info(f"File with programs paths created at '{target_programs_paths_file}'")
    logger.info("Programs installed successfully!")
