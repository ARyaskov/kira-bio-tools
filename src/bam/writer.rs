use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_sam as sam;
use noodles_sam::alignment::io::Write as AlignmentWrite;
use std::fs::File;
use std::path::Path;

pub enum OutputKind { Bam, Cram }

pub struct BamWriter {
    sink: Box<dyn BamSink + Send>,
    pub header: sam::Header,
}

trait BamSink {
    fn write_record(&mut self, header: &sam::Header, rec: &bam::Record) -> Result<()>;
    fn try_finish(&mut self, header: &sam::Header) -> Result<()>;
}

struct WrappedBam { w: Option<bam::io::Writer<noodles_bgzf::io::writer::Writer<File>>> }

impl BamSink for WrappedBam {
    fn write_record(&mut self, header: &sam::Header, rec: &bam::Record) -> Result<()> {
        self.w.as_mut().unwrap().write_record(header, rec).context("write BAM record")
    }
    fn try_finish(&mut self, _header: &sam::Header) -> Result<()> {
        if let Some(mut w) = self.w.take() { w.try_finish().context("finalize BAM")?; }
        Ok(())
    }
}

struct WrappedCram { w: Option<cram::io::Writer<File>> }

impl BamSink for WrappedCram {
    fn write_record(&mut self, header: &sam::Header, rec: &bam::Record) -> Result<()> {
        self.w.as_mut().unwrap().write_alignment_record(header, rec as &dyn sam::alignment::Record)
            .context("write CRAM record")
    }
    fn try_finish(&mut self, header: &sam::Header) -> Result<()> {
        if let Some(mut w) = self.w.take() { w.try_finish(header).context("finalize CRAM")?; }
        Ok(())
    }
}

impl BamWriter {
    pub fn create<P: AsRef<Path>>(p: P, header: &sam::Header) -> Result<Self> {
        Self::create_kind(p, header, OutputKind::Bam)
    }

    pub fn create_kind<P: AsRef<Path>>(p: P, header: &sam::Header, kind: OutputKind) -> Result<Self> {
        let f = File::create(p.as_ref()).with_context(|| format!("create {:?}", p.as_ref()))?;
        let sink: Box<dyn BamSink + Send> = match kind {
            OutputKind::Bam => {
                let mut w = bam::io::Writer::new(f);
                w.write_header(header).context("write BAM header")?;
                Box::new(WrappedBam { w: Some(w) })
            }
            OutputKind::Cram => {
                let builder = cram::io::writer::Builder::default()
                    .preserve_read_names(true);
                let mut w = builder.build_from_writer(f);
                w.write_file_definition().context("write CRAM file definition")?;
                w.write_file_header(header).context("write CRAM header")?;
                Box::new(WrappedCram { w: Some(w) })
            }
        };
        Ok(Self { sink, header: header.clone() })
    }

    pub fn write_record(&mut self, rec: &bam::Record) -> Result<()> {
        self.sink.write_record(&self.header, rec)
    }

    pub fn finish(mut self) -> Result<()> {
        self.sink.try_finish(&self.header)
    }
}
