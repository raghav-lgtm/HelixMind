# HelixMind 🧬

A full-stack genomic analysis platform that leverages advanced AI models to predict the functional impact of DNA mutations. Built with modern web technologies and cutting-edge machine learning for bioinformatics research and clinical applications.

![GitHub repo size](https://img.shields.io/github/repo-size/raghav-lgtm/HelixMind)
![GitHub last commit](https://img.shields.io/github/last-commit/raghav-lgtm/HelixMind)
![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)
![Next.js](https://img.shields.io/badge/Next.js-14+-black.svg)
![FastAPI](https://img.shields.io/badge/FastAPI-0.100+-green.svg)

## 🌟 Features

- **AI-Powered Mutation Analysis**: Utilizes EVO2 transformer models for accurate pathogenicity prediction
- **Interactive Web Interface**: Modern React-based frontend with real-time analysis
- **RESTful API Backend**: FastAPI-powered backend with async processing
- **Comprehensive Scoring**: Multiple prediction algorithms including conservation scores
- **Clinical Data Integration**: Compare predictions with known clinical databases
- **Research-Grade Tools**: Suitable for both clinical and research applications

## 🏗️ Architecture

### Backend (Python/FastAPI)

- **EVO2 Integration**: State-of-the-art transformer model for variant effect prediction
- **Scoring Engine**: Multiple algorithms for comprehensive mutation assessment
- **Async Processing**: High-performance concurrent request handling
- **Data Pipeline**: Efficient processing of genomic variants

### Frontend (Next.js/React)

- **Modern UI**: Clean, responsive interface built with Next.js 14
- **Real-time Analysis**: Interactive mutation input and instant results
- **Data Visualization**: Clear presentation of prediction scores and confidence
- **User Experience**: Intuitive workflow for genomic analysis

## 🚀 Quick Start

### Prerequisites

- Python 3.11+
- Node.js 18+
- Git

### Backend Setup

```bash
cd HelixMind-backend
python -m venv .venv
# On Windows:
.venv\Scripts\activate
# On macOS/Linux:
source .venv/bin/activate

pip install -r requirements.txt
python main.py
```

### Frontend Setup

```bash
cd HelixMind-frontend
npm install
npm run dev
```

Visit `http://localhost:3000` to access the application.

## 📊 Technical Details

### EVO2 Model Integration

- **Model Variants**: Support for 1B, 7B, and 40B parameter models
- **Context Lengths**: 8K and 1M token configurations
- **Fine-tuning**: Optimized for genomic sequence analysis
- **Performance**: GPU-accelerated inference for real-time predictions

### API Endpoints

```
POST /analyze/variant    # Single variant analysis
POST /analyze/batch      # Batch variant processing
GET  /models/info        # Available model information
GET  /health             # Service health check
```

### Scoring Metrics

- **Pathogenicity Score**: 0-1 scale indicating disease-causing potential
- **Conservation Score**: Evolutionary conservation analysis
- **Confidence Interval**: Statistical confidence in predictions
- **Clinical Correlation**: Comparison with ClinVar database

## 🔬 Use Cases

- **Clinical Diagnostics**: Support clinical decision-making for genetic variants
- **Research Analysis**: High-throughput variant screening for research studies
- **Drug Discovery**: Identify potential therapeutic targets through mutation analysis
- **Population Genomics**: Large-scale analysis of genetic variation in populations

## 🛠️ Development

### Project Structure

```
HelixMind/
├── HelixMind-backend/          # Python FastAPI backend
│   ├── evo2/                   # EVO2 model implementation
│   ├── main.py                 # API server
│   └── requirements.txt        # Python dependencies
├── HelixMind-frontend/         # Next.js frontend
│   ├── src/                    # React components
│   ├── public/                 # Static assets
│   └── package.json            # Node dependencies
└── README.md                   # Project documentation
```

### Testing

```bash
# Backend tests
cd HelixMind-backend
python -m pytest test/

# Frontend tests
cd HelixMind-frontend
npm test
```

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **EVO2 Model**: Based on the revolutionary work by [EVO2 Research Team](https://github.com/evo2/evo2)
- **FastAPI**: Modern Python web framework
- **Next.js**: React framework for production
- **Bioinformatics Community**: For advancing genomic research

## 📞 Contact

For questions or collaboration opportunities, please reach out through GitHub issues.

---

_Built with ❤️ for advancing genomic medicine and bioinformatics research_
