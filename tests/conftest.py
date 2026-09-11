"""pytest configuration: repository root on sys.path; every test run logged to outputs/logs/."""
import datetime as dt
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
FIXTURES = Path(__file__).parent / "fixtures"


def pytest_configure(config):
    logs = ROOT / "outputs" / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    config.option.log_file = str(logs / f"pytest_{stamp}.log")
    config.option.log_file_level = "INFO"
