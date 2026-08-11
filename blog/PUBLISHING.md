# Publishing the v2 CXSMILES follow-up post

All deliverables live in this `blog/` folder. Nothing here is intended to be
part of the crate's library; it is the article and its supporting assets.

## What is here

| File | Purpose |
|---|---|
| `v2_followup.md` | The post in Markdown. Use this for metabolomics.blog (or any Mermaid-rendering Markdown platform). |
| `v2_followup_wordpress.html` | A self-contained WordPress-ready HTML document. All 15 figures are embedded as base64, so this single file renders fully on its own. |
| `v2_followup_wordpress.template.html` | Source template for the above, with `{{BASE64:...}}` placeholders. Only edit this, then regenerate the final HTML with `blog/build_wordpress_html.py`. |
| `v2_followup_preview.png` | A full-page screenshot of the rendered HTML (1904x8639) for a quick visual check before publishing. |
| `v2_figures/*.png` | The 15 individual molecule depictions, for uploading to a WordPress media library if you prefer that to base64. |
| `validate_blog_examples.py` | Reproducible check that every CXSMILES string in the post parses in RDKit and CDK. |

## The two publish paths

### Markdown (metabolomics.blog)

1. Open `v2_followup.md`.
2. Upload the 15 files in `v2_figures/` to the media library.
3. Replace each `v2_figures/fig_*.png` image path with the uploaded URL (or keep the relative paths if the platform serves the folder).
4. Paste the Markdown into the post editor. The Mermaid diagram renders client-side; metabolomics.blog already does this.

### WordPress

Two options, both from `v2_followup_wordpress.html`:

- **Classic editor:** switch to the HTML tab and paste the whole file.
- **Block editor:** add a Custom HTML block and paste the whole file.

The document already carries its own `<h1>/<h2>/<h3>` headings, tables, code
blocks, the inline SVG design-rule diagram, and all 15 figures as base64, so
no separate image upload is needed. It is about 240 KB, well under normal
WordPress post size limits.

If your WordPress install strips `data:` URIs, upload the 15 PNGs in
`v2_figures/` to the media library instead and swap each
`src="data:image/png;base64,..."` for the uploaded media URL.

## Regenerating the WordPress HTML

Only the template should be edited by hand. Rebuild the final HTML after a
template change with:

```bash
python blog/build_wordpress_html.py
```

This inlines the 15 figures from `v2_figures/` and the design-rule SVG into
`v2_followup_wordpress.html`.

## Verifying the examples

The post claims every CXSMILES string parses in both CDK and RDKit. Check it
reproducibly with (requires RDKit):

```bash
uv run --python 3.12 --with rdkit python blog/validate_blog_examples.py --cdk
```

Expected output:

```
blog examples: 15
RDKit: 15/15 passed
CDK  : 15/15 passed
RESULT: PASS
```

## Style rules honored

- No em-dashes or en-dashes anywhere in the article (checked with grep).
- First-person, plain-language voice consistent with the earlier post.
- Every CXSMILES string is the literal output of the `lipid_notation`
  converter, verified byte-for-byte, not a hand-written approximation.
