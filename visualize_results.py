"""Visualize non-Gaussianity test results."""

import numpy as np
import matplotlib.pyplot as plt

# Results from the corrected test
metrics = ['Counts\nVariance', 'Counts\nSkewness', 'Euler\nCharacteristic', 'Peak\nCount']
real_values = [0.872, 0.456, 1.0, 1801]
mock_means = [1.049, 0.720, 1.0, 1744.6]
mock_stds = [0.030, 0.046, 0.0, 30.526]
anomalous = [True, True, False, False]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for idx, (metric, real, mock_mean, mock_std, is_anom) in enumerate(zip(metrics, real_values, mock_means, mock_stds, anomalous)):
    ax = axes[idx]

    # Plot mock distribution
    x_mock = [mock_mean - 3*mock_std, mock_mean + 3*mock_std] if mock_std > 0 else [mock_mean - 0.1, mock_mean + 0.1]
    ax.axvline(mock_mean, color='blue', linewidth=2, label='ΛCDM Mock Mean', linestyle='--')

    if mock_std > 0:
        ax.axvspan(mock_mean - 2*mock_std, mock_mean + 2*mock_std, alpha=0.2, color='blue', label='±2σ (ΛCDM)')
        ax.axvspan(mock_mean - mock_std, mock_mean + mock_std, alpha=0.3, color='blue')

    # Plot real value
    color = 'red' if is_anom else 'green'
    ax.axvline(real, color=color, linewidth=3, label='Real Data', linestyle='-')

    # Styling
    ax.set_title(metric, fontsize=12, fontweight='bold')
    ax.set_xlabel('Value', fontsize=10)
    ax.set_yticks([])
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='x')

    # Add text annotation
    if is_anom:
        z_score = abs(real - mock_mean) / mock_std if mock_std > 0 else 0
        ax.text(0.98, 0.95, f'⚠️ ANOMALOUS\n({z_score:.1f}σ)',
                transform=ax.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='red', alpha=0.2))
    else:
        ax.text(0.98, 0.95, '✅ Consistent',
                transform=ax.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='green', alpha=0.2))

    # Set xlim to show both values clearly
    if mock_std > 0:
        ax.set_xlim(min(real, mock_mean - 3*mock_std) * 0.95,
                    max(real, mock_mean + 3*mock_std) * 1.05)

plt.suptitle('Non-Gaussianity Test Results: Real Data vs ΛCDM Mocks\n' +
             '2/4 metrics show anomalies',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('non_gaussianity_results.png', dpi=150, bbox_inches='tight')
print("Saved: non_gaussianity_results.png")
plt.close()

# Summary figure
fig, ax = plt.subplots(figsize=(10, 6))

summary_data = {
    'All Bugs Fixed': 4,
    'Anomalous Metrics': 2,
    'Consistent Metrics': 2
}

bars = ax.bar(range(len(summary_data)), list(summary_data.values()),
              color=['green', 'orange', 'blue'], alpha=0.7, edgecolor='black', linewidth=2)

ax.set_xticks(range(len(summary_data)))
ax.set_xticklabels(list(summary_data.keys()), fontsize=11, fontweight='bold')
ax.set_ylabel('Count', fontsize=12, fontweight='bold')
ax.set_title('Implementation Status & Results Summary', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for bar in bars:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{int(height)}',
            ha='center', va='bottom', fontsize=16, fontweight='bold')

plt.tight_layout()
plt.savefig('implementation_summary.png', dpi=150, bbox_inches='tight')
print("Saved: implementation_summary.png")
plt.close()
