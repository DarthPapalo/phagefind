#!/usr/bin/env python3

############################################################################
# PhageFind - Copyright (C) 2026  Pablo (DarthPapalo) Vidal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# See file LICENSE for details.
############################################################################

import datetime
import logging
import sys
from pathlib import Path
from typing import Never, cast

import clapy
from clapy.parsed_command import ParsedCommand
from rich.logging import RichHandler

from _commons import PROGRAM_DESC, PROGRAM_NAME, VERSION, console, ensure_dir
from commands.analysis import ANALYSIS_COMMAND_NAME, AnalysisMode, run_analysis
from commands.analysis import cli as analysis_command
from commands.download_phage_data import (
    DOWNLOAD_DATA_COMMAND_NAME,
    get_data,
)
from commands.download_phage_data import (
    cli as download_phage_data_command,
)
from commands.install_programs import (
    INSTALL_PROGRAMS_COMMAND_NAME,
    run_micromamba,
)
from commands.install_programs import (
    cli as install_programs_command,
)
from commands.visualization import (
    VISUALIZATION_COMMAND_NAME,
)
from commands.visualization import (
    cli as visualization_command,
)
from commands.visualization import (
    run as run_visualization,
)
from program_config import (
    DEFAULT_CONFIG_FILE,
    DEFAULT_PROGRAMS_PATHS_FILE,
    PipelineConfig,
    ProgramsPaths,
    load_config,
    load_programs_paths,
)

# ========= CLI Setup =========
# fmt: off
cli: clapy.Command = (
    clapy.Command("phagefind.py")
        .help(f"{PROGRAM_NAME} - {PROGRAM_DESC}")
        .arguments([
                clapy.Arg("config-path")
                    .long("--config")
                    .short("-c")
                    .help("Custom configuration file to be used by the program")
                    .value_parser(Path)
                    .default(DEFAULT_CONFIG_FILE)
                    .propagate(),
                clapy.Arg("programs-paths")
                    .long("--programs-paths")
                    .short("-p")
                    .help("Path to the TOML file with the programs paths, generated automatically with install-programs command.")
                    .value_parser(Path)
                    .default(DEFAULT_PROGRAMS_PATHS_FILE)
                    .propagate(),
                clapy.Arg("debug")
                    .long("--debug")
                    .help("Enables debug level logging and other utilities.")
                    .flag()
                    .propagate(),
        ])
)
# fmt: on

# Common argument for all commands that generate data
output_argument: clapy.Arg = (
    clapy.Arg("output-path")
    .help("Path to the output directory where the results will be stored.")
    .long("--output")
    .short("-o")
    .value_parser(Path)
    .propagate()
)

cli.subcommand(download_phage_data_command.argument(output_argument))
cli.subcommand(analysis_command.argument(output_argument))
cli.subcommand(visualization_command.argument(output_argument))
cli.subcommand(install_programs_command)

version_command: clapy.Command = clapy.Command("version").help(
    "Shows the pipeline version and exits."
)
cli.subcommand(version_command)


def init_logging(level: int | str, logs_dir: Path) -> logging.Logger:
    logger: logging.Logger = logging.getLogger(PROGRAM_NAME)

    logs_dir.mkdir(parents=True, exist_ok=True)

    logger.setLevel(level)

    # ====== File logging handler ======
    file_handler = logging.FileHandler(
        logs_dir / f"{PROGRAM_NAME}_{datetime.datetime.now():%Y-%m-%d_%H-%M-%S}.log"
    )
    file_handler.setFormatter(
        logging.Formatter(
            fmt="[%(asctime)s] - %(levelname)s : %(message)s",
            datefmt="%d-%m-%Y %H:%M:%S",
        )
    )

    # ====== Console logging handler ======
    console_handler = RichHandler(console=console, show_path=False)

    # Add the different handlers
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def version() -> Never:
    """
    Prints the program version number and exits.
    """
    console.print(
        f"[b green]{PROGRAM_NAME}[/b green] [bright_black]-[/bright_black] [i bright_cyan]Version {VERSION}[/i bright_cyan]",
        highlight=False,
    )
    sys.exit(0)


def main() -> Never:
    parsed: clapy.ParsedCommand = cli.parse()

    if parsed.subcommand_name() == "version":
        version()

    elif parsed.subcommand_name() in (
        ANALYSIS_COMMAND_NAME,
        DOWNLOAD_DATA_COMMAND_NAME,
        VISUALIZATION_COMMAND_NAME,
        INSTALL_PROGRAMS_COMMAND_NAME,
    ):
        logger: logging.Logger = init_logging(
            logging.DEBUG if parsed.get_flag("debug") else logging.INFO,
            Path.cwd(),
        )

        if parsed.subcommand_name() == INSTALL_PROGRAMS_COMMAND_NAME:
            parsed_install_programs_command: clapy.ParsedCommand = parsed.subcommand()
            run_micromamba(
                logger,
                parsed_install_programs_command.get_one("programs-dir"),
                parsed_install_programs_command.get_one("programs-paths-file"),
            )

        else:
            config_path: Path = parsed.get_one("config-path")
            programs_paths: Path = parsed.get_one("programs-paths")

            loaded_config: PipelineConfig = load_config(config_path)
            loaded_programs_paths: ProgramsPaths = load_programs_paths(programs_paths)

            logger.info(f"Executing: {' '.join(sys.argv[1:])}")

            # Output dir is common to all of this subcommands
            output_dir: Path = ensure_dir(
                logger,
                cast(Path, parsed.subcommand().get_one("output-path")),
                False,
            )

            if parsed.subcommand_name() == ANALYSIS_COMMAND_NAME:
                assert parsed.subcommand().subcommand_name() in [
                    "reads",
                    "assembly",
                ]
                parsed_analysis_command: ParsedCommand = (
                    parsed.subcommand().subcommand()
                )
                run_analysis(
                    AnalysisMode.READS
                    if parsed.subcommand().subcommand_name() == "reads"
                    else AnalysisMode.ASSEMBLY,
                    logger,
                    loaded_config,
                    loaded_programs_paths,
                    output_dir,
                    parsed_analysis_command,
                )

            elif parsed.subcommand_name() == VISUALIZATION_COMMAND_NAME:
                parsed_visualization_command: clapy.ParsedCommand = parsed.subcommand()
                run_visualization(
                    logger,
                    loaded_config.visualization,
                    output_dir,
                    parsed_visualization_command.get_one("features-dir"),
                    parsed_visualization_command.get_many("features-ids", force=True),
                    parsed_visualization_command.get_many("graphs", force=True),
                )

            elif parsed.subcommand_name() == DOWNLOAD_DATA_COMMAND_NAME:
                parsed_download_phage_data_command: clapy.ParsedCommand = (
                    parsed.subcommand()
                )
                get_data(
                    logger,
                    loaded_programs_paths,
                    output_dir,
                    list(
                        parsed_download_phage_data_command.get_many(
                            "data-sources", force=True
                        )
                    ),
                    list(
                        parsed_download_phage_data_command.get_many(
                            "data-types", force=True
                        )
                    ),
                    parsed_download_phage_data_command.get_flag("keep-individual"),
                    parsed_download_phage_data_command.get_flag("verify"),
                )
            else:
                logger.critical("Couldn't parse subcommand.")
                sys.exit(1)
    else:
        # Subcommand is mandatory, help should be displayed instead
        print("Couldn't parse subcommand.")
        sys.exit(1)

    sys.exit(0)

if __name__ == "__main__":
    main()