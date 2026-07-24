import AppKit
import Foundation
import WebKit

func parseCanvasSize(from svgText: String) -> CGSize {
    let patterns = [
        #"viewBox\s*=\s*"[^"]*\s+[^"]*\s+([0-9.]+)\s+([0-9.]+)""#,
        #"width\s*=\s*"([0-9.]+)".*height\s*=\s*"([0-9.]+)""#,
    ]

    for pattern in patterns {
        if let regex = try? NSRegularExpression(pattern: pattern, options: []) {
            let range = NSRange(svgText.startIndex..<svgText.endIndex, in: svgText)
            if let match = regex.firstMatch(in: svgText, options: [], range: range),
               let wRange = Range(match.range(at: 1), in: svgText),
               let hRange = Range(match.range(at: 2), in: svgText),
               let width = Double(svgText[wRange]),
               let height = Double(svgText[hRange]) {
                return CGSize(width: width, height: height)
            }
        }
    }
    return CGSize(width: 1600, height: 1200)
}

final class AppDelegate: NSObject, NSApplicationDelegate, WKNavigationDelegate {
    let inputURL: URL
    let outputURL: URL
    let size: CGSize
    var window: NSWindow?
    var webView: WKWebView?

    init(inputURL: URL, outputURL: URL, size: CGSize) {
        self.inputURL = inputURL
        self.outputURL = outputURL
        self.size = size
    }

    func applicationDidFinishLaunching(_ notification: Notification) {
        let config = WKWebViewConfiguration()
        let frame = CGRect(origin: .zero, size: size)
        let webView = WKWebView(frame: frame, configuration: config)
        webView.navigationDelegate = self
        webView.setValue(false, forKey: "drawsBackground")
        self.webView = webView

        let window = NSWindow(
            contentRect: frame,
            styleMask: [.borderless],
            backing: .buffered,
            defer: false
        )
        window.contentView = webView
        window.orderOut(nil)
        self.window = window

        webView.loadFileURL(inputURL, allowingReadAccessTo: inputURL.deletingLastPathComponent())
    }

    func webView(_ webView: WKWebView, didFinish navigation: WKNavigation!) {
        DispatchQueue.main.asyncAfter(deadline: .now() + 0.8) {
            let config = WKSnapshotConfiguration()
            config.rect = CGRect(origin: .zero, size: self.size)
            webView.takeSnapshot(with: config) { image, error in
                if let error {
                    fputs("snapshot error: \(error)\n", stderr)
                    NSApp.terminate(nil)
                    return
                }
                guard let image,
                      let tiff = image.tiffRepresentation,
                      let rep = NSBitmapImageRep(data: tiff),
                      let png = rep.representation(using: .png, properties: [:]) else {
                    fputs("snapshot error: failed to create PNG data\n", stderr)
                    NSApp.terminate(nil)
                    return
                }
                do {
                    try png.write(to: self.outputURL)
                } catch {
                    fputs("write error: \(error)\n", stderr)
                }
                NSApp.terminate(nil)
            }
        }
    }

    func webView(_ webView: WKWebView, didFail navigation: WKNavigation!, withError error: Error) {
        fputs("navigation error: \(error)\n", stderr)
        NSApp.terminate(nil)
    }

    func webView(_ webView: WKWebView, didFailProvisionalNavigation navigation: WKNavigation!, withError error: Error) {
        fputs("provisional navigation error: \(error)\n", stderr)
        NSApp.terminate(nil)
    }
}

@main
struct RenderSVGToPNG {
    static func main() throws {
        guard CommandLine.arguments.count >= 3 else {
            fputs("usage: render_svg_to_png.swift input.svg output.png\n", stderr)
            exit(2)
        }

        let inputURL = URL(fileURLWithPath: CommandLine.arguments[1])
        let outputURL = URL(fileURLWithPath: CommandLine.arguments[2])
        let svgText = try String(contentsOf: inputURL, encoding: .utf8)
        let size = parseCanvasSize(from: svgText)

        let app = NSApplication.shared
        app.setActivationPolicy(.prohibited)
        let delegate = AppDelegate(inputURL: inputURL, outputURL: outputURL, size: size)
        app.delegate = delegate
        app.run()
    }
}
