from __future__ import annotations

from pathlib import Path
import tempfile

import streamlit as st

from web_runner import run_pipeline_for_vcf

st.set_page_config(page_title="PRS Pipeline Web", page_icon="🧬", layout="centered")

st.title("🧬 Prostate Cancer PRS Pipeline")
st.write("Sube un archivo `vcf.gz`, ejecuta el pipeline y descarga resultados.")

uploaded = st.file_uploader("Archivo VCF comprimido", type=["gz"], help="Debe ser un archivo .vcf.gz")

if uploaded is not None and not uploaded.name.endswith(".vcf.gz"):
    st.error("El nombre del archivo debe terminar en .vcf.gz")

run_button = st.button("Ejecutar pipeline", type="primary", disabled=uploaded is None)

if run_button and uploaded is not None:
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = Path(tmpdir) / uploaded.name
        input_path.write_bytes(uploaded.getvalue())

        with st.spinner("Ejecutando pipeline... esto puede tardar varios minutos"):
            try:
                result = run_pipeline_for_vcf(input_path)
            except Exception as exc:
                st.error(f"Error al ejecutar pipeline: {exc}")
            else:
                st.success(f"Pipeline completado. Job ID: {result.job_id}")
                log_text = result.log_path.read_text(encoding="utf-8")
                st.subheader("Logs")
                st.text_area("Salida del pipeline", value=log_text, height=320)

                st.download_button(
                    "Descargar log",
                    data=log_text,
                    file_name=f"{result.job_id}.log",
                    mime="text/plain",
                )

                if result.outputs:
                    st.subheader("Carpetas de salida")
                    for output in result.outputs:
                        st.code(str(output))
                else:
                    st.info("No se encontraron carpetas de salida esperadas (results_prs / 05_final).")
