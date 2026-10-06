To: Marianne Mueller <marianne.mueller@math.ethz.ch>
Subject: iarm 0.4.3: partgam_DIF() and partgam_LD() adjusted p-values are Bonferroni whatever p.adj is

Dear Marianne,

I maintain easyRasch2, which uses iarm for partial gamma and item-restscore
statistics, and I think I have found a problem with the adjusted p-values in
partgam_DIF() and partgam_LD() (iarm 0.4.3).

**Symptom.** The adjusted p-value column (labelled padj.BH by default) is a
Bonferroni correction, min(1, n × p), whatever method is passed in `p.adj`.
In the attached example with six items, the "BH" column equals k × p for
every item, and `p.adj = "holm"` returns an identical column. A real
Benjamini-Hochberg adjustment of the same p-values gives values about half
as large for the smallest p-values.

**Cause.** Both functions adjust each p-value on its own inside the item (or
item-pair) loop:

    pkorr <- p.adjust(pvalue, method = padj, n = l * k)            # partgam_DIF
    pkorr <- p.adjust(pvalue, method = padj, n = (k * (k - 1)))    # partgam_LD

With a single p-value, the step-up and step-down methods have nothing to rank
against, so BH, Holm, Hochberg and Bonferroni all return min(1, n × p). The
significance symbols are derived from this value too. item_restscore(),
out_infit() and boot_fit() adjust the whole vector at once and give the
intended result.

**Possible fix.** Collect the unadjusted p-values in the loop and call
p.adjust() once afterwards: over the items (times the number of exogenous
variables) in partgam_DIF(), and over all k(k − 1) tests in both tables
together in partgam_LD(), which is the family the current `n` describes. The
symbols would then be computed from the adjusted values.

In easyRasch2 I now apply the BH adjustment to iarm's unadjusted p-values
myself, so nothing is urgent on my side, but other users will be reading the
column as BH.

The attached reprex.R reproduces the behaviour with iarm alone.

Thank you for iarm, which easyRasch2 relies on for several of its analyses.

Best regards,
Magnus Johansson
Karolinska Institutet

Attachment: reprex.R
