"""Load the separately pinned standard-library trace arithmetic."""
from pathlib import Path
import importlib.util, os
path=Path(os.environ.get('SO11_REVIEW_PATH',str(Path(__file__).resolve().parent/'inputs/so11_derivative_review.py')))
spec=importlib.util.spec_from_file_location('pinned_independent_trace_ops',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
QUERIES=module.QUERIES
canonical=module.canonical
project=module.project
action=module.action
