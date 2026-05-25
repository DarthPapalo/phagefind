from datetime import datetime
from pathlib import Path
from typing import Literal

from jinja2 import Environment

from ._commons import AnalysisResults

TEMPLATE_LIGHT = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Analysis Report</title>
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
        body { font-family: 'Inter', Arial, sans-serif; margin: 40px; background: linear-gradient(135deg, #f5f7fa 0%, #e4e8ec 100%); min-height: 100vh; }
        .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px 40px; border-radius: 16px; margin-bottom: 30px; }
        .header h1 { color: white; margin: 0; font-size: 2.2em; font-weight: 700; }
        .header .subtitle { font-size: 0.9em; color: #e2e8f0; margin-left: 15px; font-weight: 400; }
        h2 { color: #2d3748; border-bottom: 2px solid #e2e8f0; padding-bottom: 12px; margin-top: 30px; font-weight: 600; }
        .info-section { background: white; padding: 25px; border-radius: 12px; margin-bottom: 20px; box-shadow: 0 4px 15px rgba(0,0,0,0.05); }
        .info-grid { display: grid; grid-template-columns: 120px 1fr; gap: 12px; }
        .info-label { font-weight: 600; color: #4a5568; }
        .info-value { color: #2d3748; }
        .badge { display: inline-block; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 4px 12px; border-radius: 20px; font-size: 0.85em; font-weight: 500; }
        .hit-section { background: #f8fafc; margin: 20px 0; padding: 0; border-radius: 16px; box-shadow: 0 4px 20px rgba(0,0,0,0.08); overflow: hidden; border: 1px solid #e2e8f0; }
        .feature-header { background: linear-gradient(135deg, #d7dcf3 0%, #f8fafc 100%); padding: 20px 25px; border-bottom: 1px solid #e2e8f0; }
        .feature-name { display: flex; font-size: 1.4em; font-weight: 700; color: #2d3748; margin: 0; }
        .feature-name a { color: #667eea; text-decoration: none; }
        .feature-name a:hover { color: #764ba2; text-decoration: underline; }
        .hit-details { padding: 20px 25px; }
        .hit-details-title { font-size: 1em; font-weight: 600; color: #667eea; margin-bottom: 15px; text-transform: uppercase; letter-spacing: 0.5px; }
        .detail-row { display: flex; padding: 8px 0; border-bottom: 1px solid #f0f0f0; }
        .detail-row:last-child { border-bottom: none; }
        .detail-label { width: 180px; font-weight: 600; color: #718096; flex-shrink: 0; }
        .detail-value { color: #2d3748; }
        .detail-row.phage-targets { background: linear-gradient(135deg, #c6f6d5 0%, #9ae6b4 100%); border-radius: 8px; padding: 12px 15px; margin-top: 10px; border: none; }
        .detail-row.phage-targets .detail-label { color: #22543d; width: auto; flex-shrink: 0; margin-right: 10px; font-weight: 700; }
        .detail-row.phage-targets .detail-value { color: #22543d; font-weight: 500; }
        .unsuccessful-tag { font-size: 0.7em; color: #e53e3e; margin-left: 12px; background: #fed7d7; padding: 3px 9px; border-radius: 12px; font-weight: 500; }
        .features-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
        .features-grid .hit-section { margin: 0; }
        @media (max-width: 900px) { .features-grid { grid-template-columns: 1fr; } }
    </style>
</head>
<body>
    <div class="header">
        <h1>Analysis Report <span class="subtitle">{{ report_date }}</span></h1>
    </div>

    <h2>General Information</h2>
    <div class="info-section">
        <div class="info-grid">
            <span class="info-label">Mode:</span>
            <span class="info-value"><span class="badge">{% if mode.value == 0 %}Reads analysis{% elif mode.value == 1 %}Assembly analysis{% endif %}</span></span>
            <span class="info-label">Input files:</span>
            <span class="info-value">{{ input_files | map('string') | join(', ') }}</span>
        </div>
    </div>

    <h2>Phages detected</h2>
    <div class="features-grid">
    {% for feature_id, hits in hits_per_feature.items() | sort %}
    {% set best_hit = hits | selectattr('evalue') | sort(attribute='evalue') | first %}
    {% if not best_hit %}{% set best_hit = hits | sort(attribute='bitscore', reverse=True) | first %}{% endif %}
    <div class="hit-section">
        <div class="feature-header">
            <div class="feature-name">
            {% if database == 'ENA' %}
            <a href="https://www.ebi.ac.uk/ena/browser/view/{{ feature_id }}" target="_blank">{{ feature_id }}</a>
            {% elif database == 'NCBI' %}
            <a href="https://www.ncbi.nlm.nih.gov/nuccore/{{ feature_id }}/" target="_blank">{{ feature_id }}</a>
            {% else %}
            {{ feature_id }}
            {% endif %}
            {% if feature_id not in succesful_differences_analysis_features %}<span class="unsuccessful-tag">Differences analysis failed</span>{% endif %}
            </div>
        </div>
        <div class="hit-details">
            <div class="hit-details-title">Best Hit</div>
            <div class="detail-row"><span class="detail-label">Genome Region:</span><span class="detail-value">{{ best_hit.genome_region }}</span></div>
            <div class="detail-row"><span class="detail-label">Feature Location:</span><span class="detail-value">{{ best_hit.feature_loc }}</span></div>
            <div class="detail-row"><span class="detail-label">Genome Region Location:</span><span class="detail-value">{{ best_hit.genome_region_loc }}</span></div>
            <div class="detail-row"><span class="detail-label">Bitscore:</span><span class="detail-value">{{ best_hit.bitscore }}</span></div>
            <div class="detail-row"><span class="detail-label">E-value:</span><span class="detail-value">{{ best_hit.evalue }}</span></div>
            {% if phage_targets.get(feature_id) %}
            <div class="detail-row phage-targets"><span class="detail-label">Phage targets:</span><span class="detail-value">{{ phage_targets.get(feature_id) | join(', ') }}</span></div>
            {% endif %}
        </div>
    </div>
    {% endfor %}
    </div>
</body>
</html>
"""

TEMPLATE_DARK = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Analysis Report</title>
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
        body { font-family: 'Inter', Arial, sans-serif; margin: 40px; background: linear-gradient(135deg, #1a202c 0%, #2d3748 100%); min-height: 100vh; color: #e2e8f0; }
        .header { background: linear-gradient(135deg, #0f766e 0%, #134e4a 100%); color: white; padding: 30px 40px; border-radius: 16px; margin-bottom: 30px; }
        .header h1 { color: white; margin: 0; font-size: 2.2em; font-weight: 700; }
        .header .subtitle { font-size: 0.9em; color: #a0aec0; margin-left: 15px; font-weight: 400; }
        h2 { color: #e2e8f0; border-bottom: 2px solid #4a5568; padding-bottom: 12px; margin-top: 30px; font-weight: 600; }
        .info-section { background: #2d3748; padding: 25px; border-radius: 12px; margin-bottom: 20px; box-shadow: 0 4px 15px rgba(0,0,0,0.05); }
        .info-grid { display: grid; grid-template-columns: 120px 1fr; gap: 12px; }
        .info-label { font-weight: 600; color: #a0aec0; }
        .info-value { color: #e2e8f0; }
        .badge { display: inline-block; background: linear-gradient(135deg, #0f766e 0%, #134e4a 100%); color: white; padding: 4px 12px; border-radius: 20px; font-size: 0.85em; font-weight: 500; }
        .hit-section { background: #2d3748; margin: 20px 0; padding: 0; border-radius: 16px; box-shadow: 0 4px 20px rgba(0,0,0,0.08); overflow: hidden; border: 1px solid #4a5568; }
        .feature-header { background: linear-gradient(135deg, #1a202c 0%, #2d3748 100%); padding: 20px 25px; border-bottom: 1px solid #4a5568; }
        .feature-name { display: flex; font-size: 1.4em; font-weight: 700; color: #e2e8f0; margin: 0; }
        .feature-name a { color: #63b3ed; text-decoration: none; }
        .feature-name a:hover { color: #90cdf4; text-decoration: underline; }
        .hit-details { padding: 20px 25px; }
        .hit-details-title { font-size: 1em; font-weight: 600; color: #63b3ed; margin-bottom: 15px; text-transform: uppercase; letter-spacing: 0.5px; }
        .detail-row { display: flex; padding: 8px 0; border-bottom: 1px solid #4a5568; }
        .detail-row:last-child { border-bottom: none; }
        .detail-label { width: 180px; font-weight: 600; color: #a0aec0; flex-shrink: 0; }
        .detail-value { color: #e2e8f0; }
        .detail-row.phage-targets { background: linear-gradient(135deg, #22543d 0%, #276749 100%); border-radius: 8px; padding: 12px 15px; margin-top: 10px; border: none; }
        .detail-row.phage-targets .detail-label { color: #9ae6b4; width: auto; flex-shrink: 0; margin-right: 10px; font-weight: 700; }
        .detail-row.phage-targets .detail-value { color: #9ae6b4; font-weight: 500; }
        .unsuccessful-tag { font-size: 0.7em; color: #fc8181; margin-left: 12px; background: #742a2a; padding: 3px 9px; border-radius: 12px; font-weight: 500; }
        .features-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
        .features-grid .hit-section { margin: 0; }
        @media (max-width: 900px) { .features-grid { grid-template-columns: 1fr; } }
    </style>
</head>
<body>
    <div class="header">
        <h1>Analysis Report <span class="subtitle">{{ report_date }}</span></h1>
    </div>

    <h2>General Information</h2>
    <div class="info-section">
        <div class="info-grid">
            <span class="info-label">Mode:</span>
            <span class="info-value"><span class="badge">{% if mode.value == 0 %}Reads analysis{% elif mode.value == 1 %}Assembly analysis{% endif %}</span></span>
            <span class="info-label">Input files:</span>
            <span class="info-value">{{ input_files | map('string') | join(', ') }}</span>
        </div>
    </div>

    <h2>Phages detected</h2>
    <div class="features-grid">
    {% for feature_id, hits in hits_per_feature.items() | sort %}
    {% set best_hit = hits | selectattr('evalue') | sort(attribute='evalue') | first %}
    {% if not best_hit %}{% set best_hit = hits | sort(attribute='bitscore', reverse=True) | first %}{% endif %}
    <div class="hit-section">
        <div class="feature-header">
            <div class="feature-name">
            {% if database == 'ENA' %}
            <a href="https://www.ebi.ac.uk/ena/browser/view/{{ feature_id }}" target="_blank">{{ feature_id }}</a>
            {% elif database == 'NCBI' %}
            <a href="https://www.ncbi.nlm.nih.gov/nuccore/{{ feature_id }}/" target="_blank">{{ feature_id }}</a>
            {% else %}
            {{ feature_id }}
            {% endif %}
            {% if feature_id not in succesful_differences_analysis_features %}<span class="unsuccessful-tag">Differences analysis failed</span>{% endif %}
            </div>
        </div>
        <div class="hit-details">
            <div class="hit-details-title">Best Hit</div>
            <div class="detail-row"><span class="detail-label">Genome Region:</span><span class="detail-value">{{ best_hit.genome_region }}</span></div>
            <div class="detail-row"><span class="detail-label">Feature Location:</span><span class="detail-value">{{ best_hit.feature_loc }}</span></div>
            <div class="detail-row"><span class="detail-label">Genome Region Location:</span><span class="detail-value">{{ best_hit.genome_region_loc }}</span></div>
            <div class="detail-row"><span class="detail-label">Bitscore:</span><span class="detail-value">{{ best_hit.bitscore }}</span></div>
            <div class="detail-row"><span class="detail-label">E-value:</span><span class="detail-value">{{ best_hit.evalue }}</span></div>
            {% if phage_targets.get(feature_id) %}
            <div class="detail-row phage-targets"><span class="detail-label">Phage targets:</span><span class="detail-value">{{ phage_targets.get(feature_id) | join(', ') }}</span></div>
            {% endif %}
        </div>
    </div>
    {% endfor %}
    </div>
</body>
</html>
"""


def generate_report(
    data: AnalysisResults,
    theme: Literal["light", "dark"],
    succesful_differences_analysis_features: list[str],
    database: Literal["ENA", "NCBI"],
    output_dir: Path,
) -> None:
    env = Environment()
    template = env.from_string(TEMPLATE_LIGHT if theme == "light" else TEMPLATE_DARK)
    report_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html = template.render(
        database=database,
        report_date=report_date,
        succesful_differences_analysis_features=succesful_differences_analysis_features,
        **data,
    )
    (output_dir / "analysis_report.html").write_text(html)
