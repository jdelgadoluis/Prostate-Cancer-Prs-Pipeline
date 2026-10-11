"""
Single source of truth for the PRS risk levels, based on the population percentile.

Used by 4.0-CalcPRSypercentile.py (risk column of the summary) and heatmap.py
(colour scale and categories), so that both always agree.

    percentile < 16          -> Low
    16 <= percentile <= 84   -> Average
    84 < percentile <= 97    -> High
    percentile > 97          -> Very High
"""
import math

LOW_BELOW = 16
HIGH_ABOVE = 84
VERY_HIGH_ABOVE = 97

RISK_LEVELS = ("Low", "Average", "High", "Very High")


def classify_risk(percentile):
    """Return the risk level for a percentile (0-100), or 'N/A' if missing."""
    try:
        value = float(percentile)
    except (TypeError, ValueError):
        return "N/A"
    if math.isnan(value):
        return "N/A"
    if value > VERY_HIGH_ABOVE:
        return "Very High"
    if value > HIGH_ABOVE:
        return "High"
    if value < LOW_BELOW:
        return "Low"
    return "Average"
