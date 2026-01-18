// Make code cells collapsible and add global hide/show code button
document.addEventListener('DOMContentLoaded', function() {
    const inputWrappers = document.querySelectorAll('.jp-Cell-inputWrapper');
    console.log('Found ' + inputWrappers.length + ' input wrappers');
    
    // Create and insert the hide/show code button at the top
    const button = document.createElement('button');
    button.id = 'toggle-code-button';
    button.textContent = 'Hide Code';
    button.style.cssText = `
        position: fixed;
        top: 20px;
        left: 20px;
        z-index: 1000;
        padding: 10px 20px;
        background-color: #808080;
        color: white;
        border: none;
        border-radius: 5px;
        cursor: pointer;
        font-size: 14px;
        font-weight: 600;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
        transition: background-color 0.2s ease;
    `;
    
    button.addEventListener('mouseover', function() {
        button.style.backgroundColor = '#606060';
    });
    
    button.addEventListener('mouseout', function() {
        button.style.backgroundColor = '#808080';
    });
    
    document.body.appendChild(button);
    
    // Create and insert the table of contents button
    const tocButton = document.createElement('button');
    tocButton.id = 'toggle-toc-button';
    tocButton.textContent = 'Contents';
    tocButton.style.cssText = `
        position: fixed;
        top: 60px;
        left: 20px;
        z-index: 1000;
        padding: 10px 20px;
        background-color: #808080;
        color: white;
        border: none;
        border-radius: 5px;
        cursor: pointer;
        font-size: 14px;
        font-weight: 600;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
        transition: background-color 0.2s ease;
    `;
    
    tocButton.addEventListener('mouseover', function() {
        tocButton.style.backgroundColor = '#606060';
    });
    
    tocButton.addEventListener('mouseout', function() {
        tocButton.style.backgroundColor = '#808080';
    });
    
    document.body.appendChild(tocButton);
    
    // Create the table of contents container
    const tocContainer = document.createElement('div');
    tocContainer.id = 'toc-container';
    tocContainer.style.cssText = `
        position: fixed;
        top: 110px;
        left: 20px;
        z-index: 999;
        background-color: white;
        border: 2px solid #808080;
        border-radius: 5px;
        padding: 15px;
        max-width: 250px;
        max-height: 500px;
        overflow-y: auto;
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        display: none;
    `;
    document.body.appendChild(tocContainer);
    
    // Build table of contents from headings
    function buildTableOfContents() {
        const headings = document.querySelectorAll('h1, h2, h3, h4, h5, h6');
        const toc = [];
        
        headings.forEach((heading, index) => {
            // Skip if heading doesn't have an id, create one
            if (!heading.id) {
                heading.id = `heading-${index}`;
            }
            
            // Get text content excluding the anchor-link element
            let text = heading.textContent;
            const anchorLink = heading.querySelector('a.anchor-link');
            if (anchorLink) {
                text = text.replace(anchorLink.textContent, '').trim();
            }
            
            const level = parseInt(heading.tagName[1]);
            toc.push({
                level: level,
                text: text,
                id: heading.id
            });
        });
        
        return toc;
    }
    
    // Render table of contents as a nested tree
    function renderTableOfContents(toc) {
        tocContainer.innerHTML = '';
        
        if (toc.length === 0) {
            tocContainer.innerHTML = '<p style="color: #666; font-size: 12px;">No headings found</p>';
            return;
        }
        
        const ul = document.createElement('ul');
        ul.style.cssText = `
            list-style: none;
            padding: 0;
            margin: 0;
        `;
        
        const levelStack = [{ level: 0, ul: ul }];
        
        toc.forEach(item => {
            // Pop stack back to find parent level (level should be < item.level)
            while (levelStack.length > 1 && levelStack[levelStack.length - 1].level >= item.level) {
                levelStack.pop();
            }
            
            const parentLevel = levelStack[levelStack.length - 1];
            let targetUl = parentLevel.ul;
            let currentStackLevel = parentLevel.level;
            
            // Create intermediate levels if needed
            for (let i = currentStackLevel + 1; i < item.level; i++) {
                const intermediateUl = document.createElement('ul');
                intermediateUl.style.cssText = `
                    list-style: none;
                    padding-left: 15px;
                    margin: 0;
                `;
                
                const parentLi = document.createElement('li');
                parentLi.appendChild(intermediateUl);
                targetUl.appendChild(parentLi);
                
                levelStack.push({ level: i, ul: intermediateUl });
                targetUl = intermediateUl;
            }
            
            // Create the level list if we need to go deeper
            if (levelStack[levelStack.length - 1].level < item.level) {
                const itemUl = document.createElement('ul');
                itemUl.style.cssText = `
                    list-style: none;
                    padding-left: 15px;
                    margin: 0;
                `;
                
                const parentLi = document.createElement('li');
                parentLi.appendChild(itemUl);
                targetUl.appendChild(parentLi);
                
                levelStack.push({ level: item.level, ul: itemUl });
                targetUl = itemUl;
            } else {
                // We're at same level - use the existing ul
                targetUl = levelStack[levelStack.length - 1].ul;
            }
            
            const li = document.createElement('li');
            li.style.cssText = ``;
            
            const link = document.createElement('a');
            link.href = `#${item.id}`;
            link.textContent = item.text;
            link.style.cssText = `
                color: #0066cc;
                text-decoration: none;
                cursor: pointer;
            `;
            
            // Remove any paragraph symbols or extra decorations
            link.setAttribute('aria-label', item.text);
            
            link.addEventListener('click', function(e) {
                e.preventDefault();
                const heading = document.getElementById(item.id);
                if (heading) {
                    heading.scrollIntoView({ behavior: 'smooth' });
                    // Hide the TOC after clicking
                    tocContainer.style.display = 'none';
                }
            });
            
            link.addEventListener('mouseover', function() {
                link.style.textDecoration = 'underline';
                link.style.color = '#0052a3';
            });
            
            link.addEventListener('mouseout', function() {
                link.style.textDecoration = 'none';
                link.style.color = '#0066cc';
            });
            
            li.appendChild(link);
            targetUl.appendChild(li);
        });
        
        tocContainer.appendChild(ul);
    }
    
    // Toggle table of contents visibility
    tocButton.addEventListener('click', function() {
        if (tocContainer.style.display === 'none') {
            const toc = buildTableOfContents();
            renderTableOfContents(toc);
            tocContainer.style.display = 'block';
        } else {
            tocContainer.style.display = 'none';
        }
    });
    
    let allCollapsed = false;
    
    button.addEventListener('click', function() {
        allCollapsed = !allCollapsed;
        button.textContent = allCollapsed ? 'Show Code' : 'Hide Code';
        
        inputWrappers.forEach((wrapper) => {
            if (!wrapper.parentElement.classList.contains('jp-CodeCell')) {
                return;
            }
            
            if (allCollapsed) {
                wrapper.classList.add('collapsed');
                wrapper.classList.remove('expanded');
            } else {
                wrapper.classList.remove('collapsed');
                wrapper.classList.add('expanded');
            }
        });
    });
    
    inputWrappers.forEach((wrapper, index) => {
        const codeCell = wrapper.parentElement;
        
        // Only handle code cells
        if (!codeCell.classList.contains('jp-CodeCell')) {
            return;
        }
        
        const editor = wrapper.querySelector('.jp-CodeMirrorEditor');
        if (!editor) return;
        
        // Start in expanded state
        wrapper.classList.remove('collapsed');
        wrapper.classList.add('expanded');
        
        // Add click handler
        wrapper.addEventListener('click', function(e) {
            // Ignore clicks on the code editor itself
            if (e.target.closest('.jp-CodeMirrorEditor')) {
                return;
            }
            
            e.preventDefault();
            e.stopPropagation();
            
            const isCollapsed = wrapper.classList.contains('collapsed');
            
            if (isCollapsed) {
                wrapper.classList.remove('collapsed');
                wrapper.classList.add('expanded');
            } else {
                wrapper.classList.add('collapsed');
                wrapper.classList.remove('expanded');
            }
            
            return false;
        }, true);
    });
    
    // Add copy button to code blocks
    const highlightDivs = document.querySelectorAll('.highlight');
    highlightDivs.forEach((highlightDiv) => {
        // Create wrapper for copy button
        const wrapper = document.createElement('div');
        wrapper.className = 'highlight-wrapper';
        highlightDiv.parentNode.insertBefore(wrapper, highlightDiv);
        wrapper.appendChild(highlightDiv);
        
        // Create copy button
        const copyButton = document.createElement('button');
        copyButton.className = 'copy-button';
        copyButton.textContent = 'Copy';
        copyButton.setAttribute('title', 'Copy code to clipboard');
        
        wrapper.insertBefore(copyButton, highlightDiv);
        
        // Copy functionality
        copyButton.addEventListener('click', function() {
            const pre = highlightDiv.querySelector('pre');
            if (pre) {
                const text = pre.innerText;
                navigator.clipboard.writeText(text).then(() => {
                    // Show feedback
                    const originalText = copyButton.textContent;
                    copyButton.textContent = 'Copied!';
                    copyButton.classList.add('copied');
                    
                    setTimeout(() => {
                        copyButton.textContent = originalText;
                        copyButton.classList.remove('copied');
                    }, 2000);
                }).catch(err => {
                    console.error('Failed to copy:', err);
                });
            }
        });
    });
});


