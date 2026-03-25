"""Common utilities for UDT replication scripts."""
import json
import os
import sys
import time
import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO_ROOT, 'data', 'generated')
FIG_DIR = os.path.join(REPO_ROOT, 'manuscript', 'figures')


def save_results(filename, results, subdir=None):
    """Save results dict to JSON in data/generated/.

    Parameters
    ----------
    filename : str -- output filename (e.g., '01_vacuum_profile.json')
    results : dict -- must be JSON-serializable
    subdir : str, optional -- subdirectory within data/generated/
    """
    outdir = os.path.join(DATA_DIR, subdir) if subdir else DATA_DIR
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, filename)

    # Add metadata
    results['_metadata'] = {
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()),
        'script': os.path.basename(sys.argv[0]) if sys.argv else 'unknown',
    }

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, (np.ndarray,)):
                return obj.tolist()
            if isinstance(obj, (np.bool_,)):
                return bool(obj)
            return super().default(obj)

    with open(path, 'w') as f:
        json.dump(results, f, indent=2, cls=NumpyEncoder)
    print(f"  -> Saved: {path}")
    return path


def load_results(filename, subdir=None):
    """Load results from JSON in data/generated/."""
    outdir = os.path.join(DATA_DIR, subdir) if subdir else DATA_DIR
    path = os.path.join(outdir, filename)
    with open(path) as f:
        return json.load(f)


def pct_error(predicted, observed):
    """Percentage error: (pred - obs) / obs * 100."""
    return (predicted - observed) / observed * 100


def sigma_tension(predicted, observed, sigma):
    """Tension in units of sigma."""
    return abs(predicted - observed) / sigma


def gate_check(name, value, threshold, description=""):
    """Check a consistency gate and print result."""
    passed = abs(value) < threshold
    status = "PASS" if passed else "FAIL"
    print(f"  Gate {name}: {status} (|{value:.2e}| {'<' if passed else '>='} {threshold:.0e}) {description}")
    return passed


def format_comparison(name, predicted, observed, unit=""):
    """Format a prediction vs observation comparison."""
    err = pct_error(predicted, observed)
    u = f" {unit}" if unit else ""
    print(f"  {name}: {predicted:.6f}{u} vs {observed:.6f}{u} ({err:+.3f}%)")
    return err


def savefig(fig, name, dpi=300):
    """Save figure to manuscript/figures/."""
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ['pdf', 'png']:
        path = os.path.join(FIG_DIR, f"{name}.{ext}")
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
    print(f"  -> Saved: {os.path.join(FIG_DIR, name)}.{{pdf,png}}")
