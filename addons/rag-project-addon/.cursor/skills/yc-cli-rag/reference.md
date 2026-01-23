# yc-cli RAG API quick reference

## Endpoints

- `GET /health` -> `{"status":"ok","port":8733}`
- `GET /api/status` -> `{"count": <int>, "last_updated": <str>, ...}`
- `POST /api/index` -> `{"count": <int>, "files": [...] }`
- `GET /api/search?q=...&top_k=...` -> `[{"id","score","title","file_path"}, ...]`
- `POST /api/record` -> `{"success": true, "result": {"file_path": "...", ...}}`

## Record payload example

```json
{
  "task": "your-project",
  "title": "Decision title",
  "type": "decision",
  "conclusion": "One-line conclusion",
  "background": "Background information",
  "key_points": ["Point 1", "Point 2"],
  "tags": ""
}
```
