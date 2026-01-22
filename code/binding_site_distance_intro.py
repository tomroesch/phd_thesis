# --- Histogram parameters ---
bins = np.arange(-300, 305, 10)

# --- Activators ---
plt.figure()
plt.hist(
    tfb2.loc[is_activator, "dist_bp_signed_oriented"],
    bins=bins
)
plt.xlim(-300, 300)
plt.xlabel("TFBS center − nearest TSS (bp, strand-oriented)")
plt.ylabel("Count")
plt.title("Activator binding sites relative to TSS")
plt.tight_layout()
plt.savefig("tfbs_to_tss_activators_hist.png", dpi=200)
print("Saved tfbs_to_tss_activators_hist.png")

# --- Repressors ---
plt.figure()
plt.hist(
    tfb2.loc[is_repressor, "dist_bp_signed_oriented"],
    bins=bins
)
plt.xlim(-300, 300)
plt.xlabel("TFBS center − nearest TSS (bp, strand-oriented)")
plt.ylabel("Count")
plt.title("Repressor binding sites relative to TSS")
plt.tight_layout()
plt.savefig("tfbs_to_tss_repressors_hist.png", dpi=200)
print("Saved tfbs_to_tss_repressors_hist.png")
